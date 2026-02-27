"""RcalPySCF"""
import functools
from pathlib import Path
from typing import List, Literal, Tuple

from pyscf import dft, gto, lib, scf
from pyscf.geomopt import geometric_solver

from mcal.calculations.rcal import Rcal


print = functools.partial(print, flush=True)


class RcalPySCF(Rcal):
    """Calculate reorganization energy using PySCF."""
    HARTREE_TO_EV = 27.2114

    def __init__(
        self,
        xyz_file: str,
        osc_type: Literal['p', 'n'] = 'p',
        method: str = 'B3LYP/6-31G(d,p)',
        use_gpu: bool = False,
        ncore: int = 4,
        max_memory_gb: int = 16,
        cart: bool = False,
    ) -> None:
        super().__init__(xyz_file, osc_type, fmt='xyz')
        self._method = method
        self._use_gpu = use_gpu
        self._ncore = ncore
        self._max_memory_gb = max_memory_gb
        self._cart = cart

    def calc_reorganization(
        self,
        gau_com: str = 'g16',
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_specified_cal: List[Literal['opt_neutral', 'opt_ion', 'neutral', 'ion']] = [],
    ) -> float:
        """Calculate reorganization energy.

        Parameters
        ----------
        gau_com : str
            Ignored. Kept for interface compatibility with Rcal.
        only_read : bool
            If True, load results from existing checkpoint files without running calculations.
        is_output_detail : bool
            If True, print energy details for each step.
        skip_specified_cal : List[Literal['opt_neutral', 'opt_ion', 'neutral', 'ion']]
            Steps to skip (load from existing files instead of running).

        Returns
        -------
        float
            Reorganization energy [eV].
        """
        functional, basis = self._method.split('/', 1)

        file_path = Path(self.input_file)
        filename = file_path.stem.replace('_opt_n', '')
        directory = file_path.parent
        basename = f'{directory}/{filename}'

        atoms = self._read_xyz(self.input_file)
        energy = []

        # Step 1: opt_neutral
        mol_n = self._build_mol(atoms, charge=0, spin=0, basis=basis)
        if only_read or 'opt_neutral' in skip_specified_cal:
            pass
        else:
            print(f'running geometric optimization for {basename}_opt_n.chk')
        mol_opt_n, e0 = self._run_opt(
            mol_n, label=f'{basename}_opt_n',
            functional=functional, basis=basis,
            only_read=only_read, is_output_detail=is_output_detail,
            skip_cal='opt_neutral' in skip_specified_cal,
        )
        energy.append(e0)

        ion_charge = +1 if self.ion == 'c' else -1
        atoms_opt_n = [
            (mol_opt_n.atom_symbol(i), tuple(mol_opt_n.atom_coords(unit='Angstrom')[i]))
            for i in range(mol_opt_n.natm)
        ]

        # Step 2: ion SP at neutral geometry
        mol_ion_n = self._build_mol(atoms_opt_n, charge=ion_charge, spin=1, basis=basis)
        if only_read or 'ion' in skip_specified_cal:
            pass
        else:
            print(f'running SCF calculation for {basename}_{self.ion}.chk')
        e1 = self._run_sp(
            mol_ion_n, label=f'{basename}_{self.ion}',
            functional=functional,
            only_read=only_read, is_output_detail=is_output_detail,
            skip_cal='ion' in skip_specified_cal,
        )
        energy.append(e1)

        # Step 3: opt_ion (start from neutral-optimized geometry)
        if only_read or 'opt_ion' in skip_specified_cal:
            pass
        else:
            print(f'running geometric optimization for {basename}_opt_{self.ion}.chk')
        mol_opt_ion, e2 = self._run_opt(
            mol_ion_n, label=f'{basename}_opt_{self.ion}',
            functional=functional, basis=basis,
            only_read=only_read, is_output_detail=is_output_detail,
            skip_cal='opt_ion' in skip_specified_cal,
        )
        energy.append(e2)

        atoms_opt_ion = [
            (mol_opt_ion.atom_symbol(i), tuple(mol_opt_ion.atom_coords(unit='Angstrom')[i]))
            for i in range(mol_opt_ion.natm)
        ]

        # Step 4: neutral SP at ion geometry
        mol_n_i = self._build_mol(atoms_opt_ion, charge=0, spin=0, basis=basis)
        if only_read or 'neutral' in skip_specified_cal:
            pass
        else:
            print(f'running SCF calculation for {basename}_n.chk')
        e3 = self._run_sp(
            mol_n_i, label=f'{basename}_n',
            functional=functional,
            only_read=only_read, is_output_detail=is_output_detail,
            skip_cal='neutral' in skip_specified_cal,
        )
        energy.append(e3)

        return (energy[3] - energy[2]) + (energy[1] - energy[0])

    def _build_mol(self, atoms: List[Tuple], charge: int, spin: int, basis: str) -> gto.Mole:
        """Build molecule.

        Parameters
        ----------
        atoms : List[Tuple]
            List of atoms.
        charge : int
            Charge of the molecule.
        spin : int
            Spin of the molecule.
        basis : str
            Basis set.

        Returns
        -------
        gto.Mole
            Molecule.
        """
        mol = gto.Mole()
        mol.atom = atoms
        mol.basis = basis
        mol.charge = charge
        mol.spin = spin
        mol.verbose = 0
        mol.symmetry = False
        mol.max_memory = self._max_memory_gb * 1000
        mol.cart = self._cart
        mol.build()
        return mol

    def _build_mf(self, mol: gto.Mole, functional: str):
        """Build mean field.

        Parameters
        ----------
        mol : gto.Mole
            Molecule.
        functional : str
            Functional.

        Returns
        -------
        mf : scf.Mole
            Mean field.
        """
        if self._use_gpu:
            try:
                if functional.upper() == 'HF':
                    from gpu4pyscf.scf import hf as gpu_hf, uhf as gpu_uhf
                    mf = gpu_hf.RHF(mol) if mol.spin == 0 else gpu_uhf.UHF(mol)
                else:
                    from gpu4pyscf.dft import rks as gpu_rks, uks as gpu_uks
                    mf = gpu_rks.RKS(mol) if mol.spin == 0 else gpu_uks.UKS(mol)
                    mf.xc = functional
            except ImportError:
                print("Error: gpu4pyscf is not installed.")
                exit(1)
        else:
            if functional.upper() == 'HF':
                mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
            else:
                mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
                mf.xc = functional
        return mf

    def _run_sp(
        self,
        mol: gto.Mole,
        label: str,
        functional: str,
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_cal: bool = False,
    ) -> float:
        """Run SCF calculation.

        Parameters
        ----------
        mol : gto.Mole
            Molecule.
        label : str
            Label of the calculation.
        functional : str
            Functional.
        only_read : bool
            If True, load results from existing checkpoint files without running calculations.
        is_output_detail : bool
            If True, print energy details for each step.
        skip_cal : bool
            If True, skip the calculation.

        Returns
        -------
        float
            Energy [eV].
        """
        chkfile = f'{label}.chk'
        if only_read or skip_cal:
            energy_ev = lib.chkfile.load(chkfile, 'scf/e_tot') * self.HARTREE_TO_EV
        else:
            lib.num_threads(self._ncore)
            mf = self._build_mf(mol, functional)
            mf.chkfile = chkfile
            mf.kernel()
            if not mf.converged:
                print(f'WARNING: SCF did not converge for {label}')
            energy_ev = mf.e_tot * self.HARTREE_TO_EV

        if not only_read and not skip_cal:
            lib.chkfile.save(chkfile, 'job_status/completed', True)

        if is_output_detail:
            if not only_read and not skip_cal:
                print(f'{chkfile} calculation completed.')
            elif skip_cal:
                print(f'{chkfile} calculation skipped.')
            print(f'reading {chkfile}')
            print()
            print('--------------')
            print(' Total energy ')
            print('--------------')
            print(f'{energy_ev:12.6f} eV')
            print()

        return energy_ev

    def _run_opt(
        self,
        mol: gto.Mole,
        label: str,
        functional: str,
        basis: str,
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_cal: bool = False,
    ) -> Tuple[gto.Mole, float]:
        """Run geometric optimization.

        Parameters
        ----------
        mol : gto.Mole
            Molecule.
        label : str
            Label of the calculation.
        functional : str
            Functional.
        basis : str
            Basis set.
        only_read : bool
            If True, load results from existing checkpoint files without running calculations.
        is_output_detail : bool
            If True, print energy details for each step.
        skip_cal : bool
            If True, skip the calculation.

        Returns
        -------
        Tuple[gto.Mole, float]
            Optimized molecule and energy [eV].
        """
        xyz_file = f'{label}.xyz'
        chkfile = f'{label}.chk'
        if only_read or skip_cal:
            atoms = self._read_xyz(xyz_file)
            mol_opt = self._build_mol(atoms, mol.charge, mol.spin, basis)
            energy_ev = lib.chkfile.load(chkfile, 'scf/e_tot') * self.HARTREE_TO_EV
        else:
            lib.num_threads(self._ncore)
            mf = self._build_mf(mol, functional)
            log_ini = Path(__file__).parent / 'log.ini'
            mol_opt = geometric_solver.optimize(mf, logIni=str(log_ini))
            mf_opt = self._build_mf(mol_opt, functional)
            mf_opt.chkfile = chkfile
            mf_opt.kernel()
            if not mf_opt.converged:
                print(f'WARNING: SCF did not converge for {label}')
            energy_ev = mf_opt.e_tot * self.HARTREE_TO_EV
            self._save_xyz(mol_opt, xyz_file, label)

        if not only_read and not skip_cal:
            lib.chkfile.save(chkfile, 'job_status/completed', True)

        if is_output_detail:
            if not only_read and not skip_cal:
                print(f'{xyz_file} calculation completed.')
            elif skip_cal:
                print(f'{xyz_file} calculation skipped.')
            print(f'reading {chkfile}')
            print()
            print('--------------')
            print(' Total energy ')
            print('--------------')
            print(f'{energy_ev:12.6f} eV')
            print()

        return mol_opt, energy_ev

    @staticmethod
    def _read_xyz(xyz_file: str) -> List[Tuple]:
        """Read xyz file.

        Parameters
        ----------
        xyz_file : str
            Path to the xyz file.

        Returns
        -------
        List[Tuple]
            List of atoms.
        """
        atoms = []
        with open(xyz_file, 'r', encoding='utf-8') as f:
            f.readline()  # skip atom count
            f.readline()  # skip comment
            while True:
                line = f.readline()
                if not line:
                    break
                parts = line.strip().split()
                if len(parts) == 4:
                    symbol = parts[0]
                    x, y, z = map(float, parts[1:4])
                    atoms.append((symbol, (x, y, z)))
        return atoms

    def _save_xyz(self, mol: gto.Mole, xyz_file: str, title: str) -> None:
        """Save xyz file.

        Parameters
        ----------
        mol : gto.Mole
            Molecule.
        xyz_file : str
            Path to the xyz file.
        title : str
            Title of the molecule.
        """
        coords_ang = mol.atom_coords(unit='Angstrom')
        with open(xyz_file, 'w', encoding='utf-8') as f:
            f.write(f'{mol.natm}\n{title}\n')
            for i in range(mol.natm):
                sym = mol.atom_symbol(i)
                x, y, z = coords_ang[i]
                f.write(f'{sym}  {x:.8f}  {y:.8f}  {z:.8f}\n')
