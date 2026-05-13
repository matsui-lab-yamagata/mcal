"""RcalOrca"""
import functools
import os
from pathlib import Path
from typing import List, Literal, Optional, Tuple

from opi.core import Calculator
from opi.input.simple_keywords.task import Task
from opi.input.structures.structure import Structure
from opi.output.core import Output

from mcal.calculations.rcal import Rcal


print = functools.partial(print, flush=True)


class RcalORCA(Rcal):
    """Calculate reorganization energy using ORCA via OPI."""
    HARTREE_TO_EV = 27.2114

    def __init__(
        self,
        xyz_file: str,
        osc_type: Literal['p', 'n'] = 'p',
        method: str = 'B3LYP/6-31G(d,p)',
        ncore: int = 4,
        max_memory_gb: int = 16,
        open_mpi_path: Optional[str] = None,
    ) -> None:
        super().__init__(xyz_file, osc_type, fmt='xyz')
        if '/' not in method:
            raise ValueError(f"Method '{method}' must be in 'FUNCTIONAL/BASIS' format.")
        self._method = method
        self._ncore = ncore
        self._max_memory_gb = max_memory_gb
        self._open_mpi_path = open_mpi_path

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
            If True, load results from existing output files without running calculations.
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

        ion_charge = +1 if self.ion == 'c' else -1
        energy = []

        # Step 1: opt_neutral
        label = f'{basename}_opt_n'
        input_xyz = f'{label}_input.xyz'
        atoms = self._read_xyz(self.input_file)
        self._save_xyz(atoms, input_xyz, label)
        if only_read or 'opt_neutral' in skip_specified_cal:
            pass
        else:
            print(f'running geometric optimization for {label}.out')
        atoms_opt_n, e0 = self._run_opt(
            label=label,
            xyz_file=input_xyz,
            functional=functional,
            basis=basis,
            charge=0,
            multiplicity=1,
            only_read=only_read,
            is_output_detail=is_output_detail,
            skip_cal='opt_neutral' in skip_specified_cal,
        )
        energy.append(e0)

        # Step 2: ion SP at neutral geometry
        label = f'{basename}_{self.ion}'
        input_xyz = f'{label}_input.xyz'
        self._save_xyz(atoms_opt_n, input_xyz, label)
        moinp_n = f'{basename}_opt_n.gbw'
        if only_read or 'ion' in skip_specified_cal:
            pass
        else:
            print(f'running SCF calculation for {label}.out')
        e1 = self._run_sp(
            label=label,
            xyz_file=input_xyz,
            functional=functional,
            basis=basis,
            charge=ion_charge,
            multiplicity=2,
            only_read=only_read,
            is_output_detail=is_output_detail,
            skip_cal='ion' in skip_specified_cal,
            moinp_gbw=moinp_n,
        )
        energy.append(e1)

        # Step 3: opt_ion (start from neutral-optimized geometry)
        label = f'{basename}_opt_{self.ion}'
        input_xyz = f'{label}_input.xyz'
        self._save_xyz(atoms_opt_n, input_xyz, label)
        moinp_ion = f'{basename}_{self.ion}.gbw'
        if only_read or 'opt_ion' in skip_specified_cal:
            pass
        else:
            print(f'running geometric optimization for {label}.out')
        atoms_opt_ion, e2 = self._run_opt(
            label=label,
            xyz_file=input_xyz,
            functional=functional,
            basis=basis,
            charge=ion_charge,
            multiplicity=2,
            only_read=only_read,
            is_output_detail=is_output_detail,
            skip_cal='opt_ion' in skip_specified_cal,
            moinp_gbw=moinp_ion,
        )
        energy.append(e2)

        # Step 4: neutral SP at ion geometry
        label = f'{basename}_n'
        input_xyz = f'{label}_input.xyz'
        self._save_xyz(atoms_opt_ion, input_xyz, label)
        moinp_opt_ion = f'{basename}_opt_{self.ion}.gbw'
        if only_read or 'neutral' in skip_specified_cal:
            pass
        else:
            print(f'running SCF calculation for {label}.out')
        e3 = self._run_sp(
            label=label,
            xyz_file=input_xyz,
            functional=functional,
            basis=basis,
            charge=0,
            multiplicity=1,
            only_read=only_read,
            is_output_detail=is_output_detail,
            skip_cal='neutral' in skip_specified_cal,
            moinp_gbw=moinp_opt_ion,
        )
        energy.append(e3)

        return (energy[3] - energy[2]) + (energy[1] - energy[0])

    def _build_calculator(
        self,
        label: str,
        xyz_file: str,
        charge: int,
        multiplicity: int,
        functional: str,
        basis: str,
        task: Task,
        moinp_gbw: Optional[str] = None,
    ) -> Calculator:
        """Build an ORCA Calculator for a given step.

        Parameters
        ----------
        label : str
            Full path stem (e.g. '/path/to/mol_opt_n').
        xyz_file : str
            Path to XYZ file providing the coordinates.
        charge : int
            Molecular charge.
        multiplicity : int
            Spin multiplicity.
        functional : str
            DFT functional (e.g. 'B3LYP').
        basis : str
            Basis set (e.g. '6-31G(d,p)').
        task : Task
            Task.SP or Task.OPT.
        moinp_gbw : str, optional
            Path to GBW file for initial orbital guess.

        Returns
        -------
        Calculator
        """
        label_path = Path(label)
        calc = Calculator(
            basename=str(label_path.name),
            working_dir=label_path.parent,
            version_check=False,
        )
        calc.structure = Structure.from_xyz(xyz_file, charge=charge, multiplicity=multiplicity)
        calc.input.add_simple_keywords(functional, basis, task)
        calc.input.ncores = self._ncore
        calc.input.memory = int(self._max_memory_gb * 1024 / self._ncore)
        if moinp_gbw is not None and Path(moinp_gbw).is_file():
            calc.input.moinp = moinp_gbw
        return calc

    def _run_sp(
        self,
        label: str,
        xyz_file: str,
        functional: str,
        basis: str,
        charge: int,
        multiplicity: int,
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_cal: bool = False,
        moinp_gbw: Optional[str] = None,
    ) -> float:
        """Run or read a single-point SCF calculation.

        Parameters
        ----------
        label : str
            Full path stem of the calculation.
        xyz_file : str
            Path to XYZ file for the structure.
        functional : str
            DFT functional.
        basis : str
            Basis set.
        charge : int
            Molecular charge.
        multiplicity : int
            Spin multiplicity.
        only_read : bool
            If True, parse existing output file without running.
        is_output_detail : bool
            If True, print energy details.
        skip_cal : bool
            If True, skip calculation and read existing output.
        moinp_gbw : str, optional
            Path to GBW file for initial orbital guess.

        Returns
        -------
        float
            Energy [eV].
        """
        out_file = f'{label}.out'

        if only_read or skip_cal:
            label_path = Path(label)
            output = Output(
                basename=str(label_path.name),
                working_dir=label_path.parent,
                version_check=False,
            )
            output.parse()
            energy_hartree = output.get_final_energy()
            if energy_hartree is None:
                raise ValueError(f'Could not parse energy from {out_file}')
            energy_ev = energy_hartree * self.HARTREE_TO_EV
        else:
            calc = self._build_calculator(
                label, xyz_file, charge, multiplicity,
                functional, basis, Task.SP, moinp_gbw,
            )
            _prev_opi_mpi = os.environ.get('OPI_MPI') if self._open_mpi_path is not None else None
            if self._open_mpi_path is not None:
                os.environ['OPI_MPI'] = self._open_mpi_path
            try:
                calc.write_and_run()
            finally:
                if self._open_mpi_path is not None:
                    if _prev_opi_mpi is None:
                        del os.environ['OPI_MPI']
                    else:
                        os.environ['OPI_MPI'] = _prev_opi_mpi
            output = calc.get_output()
            output.parse()
            energy_hartree = output.get_final_energy()
            if energy_hartree is None:
                raise ValueError(f'Could not parse energy from {out_file}')
            energy_ev = energy_hartree * self.HARTREE_TO_EV

        if is_output_detail:
            if not only_read and not skip_cal:
                print(f'{out_file} calculation completed.')
            elif skip_cal:
                print(f'{out_file} calculation skipped.')
            print(f'reading {out_file}')
            print()
            print('--------------')
            print(' Total energy ')
            print('--------------')
            print(f'{energy_ev:12.6f} eV')
            print()

        return energy_ev

    def _run_opt(
        self,
        label: str,
        xyz_file: str,
        functional: str,
        basis: str,
        charge: int,
        multiplicity: int,
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_cal: bool = False,
        moinp_gbw: Optional[str] = None,
    ) -> Tuple[List[Tuple], float]:
        """Run or read a geometry optimization calculation.

        Parameters
        ----------
        label : str
            Full path stem of the calculation.
        xyz_file : str
            Path to XYZ file for the initial structure.
        functional : str
            DFT functional.
        basis : str
            Basis set.
        charge : int
            Molecular charge.
        multiplicity : int
            Spin multiplicity.
        only_read : bool
            If True, read existing output and saved XYZ without running.
        is_output_detail : bool
            If True, print energy details.
        skip_cal : bool
            If True, skip calculation and read existing output.
        moinp_gbw : str, optional
            Path to GBW file for initial orbital guess.

        Returns
        -------
        Tuple[List[Tuple], float]
            Optimized atoms list and energy [eV].
        """
        out_file = f'{label}.out'

        if only_read or skip_cal:
            label_path = Path(label)
            output = Output(
                basename=str(label_path.name),
                working_dir=label_path.parent,
                version_check=False,
            )
            output.parse()
            energy_hartree = output.get_final_energy()
            if energy_hartree is None:
                raise ValueError(f'Could not parse energy from {out_file}')
            energy_ev = energy_hartree * self.HARTREE_TO_EV
            structure = output.get_structure(index=-1)
            if structure is None:
                raise ValueError(f'Could not extract optimized geometry from {out_file}')
            atoms = [
                (str(atom.element), (atom.coordinates.x, atom.coordinates.y, atom.coordinates.z))
                for atom in structure.atoms
            ]
        else:
            calc = self._build_calculator(
                label, xyz_file, charge, multiplicity,
                functional, basis, Task.OPT, moinp_gbw,
            )
            _prev_opi_mpi = os.environ.get('OPI_MPI') if self._open_mpi_path is not None else None
            if self._open_mpi_path is not None:
                os.environ['OPI_MPI'] = self._open_mpi_path
            try:
                calc.write_and_run()
            finally:
                if self._open_mpi_path is not None:
                    if _prev_opi_mpi is None:
                        del os.environ['OPI_MPI']
                    else:
                        os.environ['OPI_MPI'] = _prev_opi_mpi
            output = calc.get_output()
            output.parse()
            energy_hartree = output.get_final_energy()
            if energy_hartree is None:
                raise ValueError(f'Could not parse energy from {out_file}')
            energy_ev = energy_hartree * self.HARTREE_TO_EV

            structure = output.get_structure(index=-1)
            if structure is None:
                raise ValueError(f'Could not extract optimized geometry from {out_file}')
            atoms = [
                (str(atom.element), (atom.coordinates.x, atom.coordinates.y, atom.coordinates.z))
                for atom in structure.atoms
            ]

        if is_output_detail:
            if not only_read and not skip_cal:
                print(f'{out_file} calculation completed.')
            elif skip_cal:
                print(f'{out_file} calculation skipped.')
            print(f'reading {out_file}')
            print()
            print('--------------')
            print(' Total energy ')
            print('--------------')
            print(f'{energy_ev:12.6f} eV')
            print()

        return atoms, energy_ev

    @staticmethod
    def _read_xyz(xyz_file: str) -> List[Tuple]:
        """Read XYZ file.

        Parameters
        ----------
        xyz_file : str
            Path to the XYZ file.

        Returns
        -------
        List[Tuple]
            List of (symbol, (x, y, z)) tuples.
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

    @staticmethod
    def _save_xyz(atoms: List[Tuple], xyz_file: str, title: str) -> None:
        """Save atoms list to XYZ file.

        Parameters
        ----------
        atoms : List[Tuple]
            List of (symbol, (x, y, z)) tuples.
        xyz_file : str
            Path to output XYZ file.
        title : str
            Title line for the XYZ file.
        """
        with open(xyz_file, 'w', encoding='utf-8') as f:
            f.write(f'{len(atoms)}\n{title}\n')
            for symbol, (x, y, z) in atoms:
                f.write(f'{symbol}  {x:.8f}  {y:.8f}  {z:.8f}\n')
