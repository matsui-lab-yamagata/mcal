"""mcal"""
import argparse
import functools
import pickle
from pathlib import Path
from time import time
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from mcal.utils.cif_reader import CifReader
from mcal.utils.gaus_log_reader import check_normal_termination
from mcal.utils.gjf_maker import GjfMaker
from mcal.calculations.hopping_mobility_model import (
    diffusion_coefficient_tensor,
    diffusion_coefficient_tensor_MC,
    diffusion_coefficient_tensor_ODE,
    marcus_rate,
    mobility_tensor
)
from mcal.calculations.rcal import Rcal
from mcal.calculations.tcal import Tcal


print = functools.partial(print, flush=True)


def main():
    """Calculate mobility tensor considering anisotropy and path continuity.

    Examples
    --------
    Basic usage:
        - Calculate p-type mobility for xxx crystal\n
        $ python hop_mcal.py xxx.cif p

        - Calculate n-type mobility for xxx crystal\n
        $ python hop_mcal.py xxx.cif n

    With resource options:
        - Use 8 CPUs and 16GB memory\n
        $ python hop_mcal.py xxx.cif p -c 8 -m 16

        - Use different calculation method (default is B3LYP/6-31G(d,p))\n
        $ python hop_mcal.py xxx.cif p -M "B3LYP/6-311G(d,p)"

    High-precision calculation:
        - Calculate all transfer integrals without speedup using moment of inertia and distance between centers of weight\n
        $ python hop_mcal.py xxx.cif p --fullcal

        - Expand calculation range to 3x3x3 supercell\n
        $ python hop_mcal.py xxx.cif p --cellsize 1

        - Expand calculation range to 5x5x5 supercell to widen transfer integral calculation range\n
        $ python hop_mcal.py xxx.cif p --cellsize 2

    Resume and save results:
        - Resume from existing calculations\n
        $ python hop_mcal.py xxx.cif p --resume

        - Save results to pickle file\n
        $ python hop_mcal.py xxx.cif p --pickle

        - Read results from existing pickle file\n
        $ python hop_mcal.py xxx_result.pkl p -rp

        - Read results from existing log files without running Gaussian\n
        $ python hop_mcal.py xxx.cif p -r

    Compare calculation methods:
        - Compare results using kinetic Monte Carlo and ODE methods\n
        $ python hop_mcal.py xxx.cif p --mc --ode
    """
    # Error range for skipping calculation of transfer integrals using moment of inertia and distance between centers of weight.
    CENTER_OF_WEIGHT_ERROR = 1.0e-7
    MOMENT_OF_INERTIA_ERROR = np.array([[1.0e-3, 1.0e-3, 1.0e-3]])

    """This code is to execute hop_mcal for command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='cif file name or pickle file name if you want to use -rp option', type=str)
    parser.add_argument('osc_type', help='organic semiconductor type', type=str)
    parser.add_argument(
        '-M', '--method',
        help='calculation method used in Gaussian calculations (default is B3LYP/6-31G(d,p))',
        type=str,
        default='B3LYP/6-31G(d,p)',
    )
    parser.add_argument('-c', '--cpu', help='setting the number of cpu (default is 4)', type=int, default=4)
    parser.add_argument(
        '-m', '--mem',
        help='setting the number of memory [GB] (default is 10 GB)',
        type=int,
        default=10,
    )
    parser.add_argument('-g', '--g09', help='use Gaussian 09 (default is Gaussian 16)', action='store_true')
    parser.add_argument('-r', '--read', help='read log files without executing Gaussian', action='store_true')
    parser.add_argument(
        '-rp', '--read_pickle',
        help='read results from existing pickle file',
        action='store_true'
    )
    parser.add_argument('-p', '--pickle', help='save to pickle the result of calculation', action='store_true')
    parser.add_argument(
        '--cellsize',
        help='number of unit cells to expand in each direction around the central unit cell '
            '(Examples: 1 creates 3x3x3, 2 creates 5x5x5 supercell (default is 2))',
        type=int,
        default=2,
    )
    parser.add_argument(
        '--fullcal',
        help='do not process for speeding up using moment of inertia and distance between centers of weight',
        action='store_true',
    )
    parser.add_argument('--mc', help='use Monte Carlo method to calculate diffusion coefficient', action='store_true')
    parser.add_argument(
        '--ode',
        help='use Ordinary Differential Equation method to calculate diffusion coefficient',
        action='store_true',
    )
    parser.add_argument(
        '--resume',
        help='resume calculation',
        action='store_true',
    )
    args = parser.parse_args()

    args.osc_type = args.osc_type.lower()

    if args.g09:
        gau_com = 'g09'
    else:
        gau_com = 'g16'

    # file info
    cif_file = Path(args.file)
    directory = cif_file.parent
    filename = cif_file.stem
    cif_path_without_ext = f'{directory}/{filename}'

    print('----------------------------------------')
    print(' mcal 0.1.1 (2026/01/08) by Matsui Lab. ')
    print('----------------------------------------')

    if args.read_pickle:
        read_pickle(args.file)
        exit()

    print(f'\nCalculate as {args.osc_type}-type organic semiconductor.')
    print(f'\nInput File Name: {args.file}')
    Tcal.print_timestamp()
    print()
    start_time = time()

    ##### Calculate reorganization energy #####
    cif_reader = CifReader(cif_path=cif_file)
    print(f'Export {cif_path_without_ext}_unit_cell.mol')
    cif_reader.export_unit_cell_file(f'{cif_path_without_ext}_unit_cell.mol', format='mol')
    print('Please verify that the created unit cell is correct.\n')
    symbols = cif_reader.unique_symbols[0]
    coordinates = cif_reader.unique_coords[0]
    coordinates = cif_reader.convert_frac_to_cart(coordinates)

    if not args.read:
        print('Create gjf for reorganization energy.')
        create_reorg_gjf(
            symbols,
            coordinates,
            filename,
            directory,
            args.cpu,
            args.mem,
            args.method,
        )

    if args.osc_type == 'p':
        rcal = Rcal(gjf_file=f'{cif_path_without_ext}_opt_n.gjf')
    elif args.osc_type == 'n':
        rcal = Rcal(gjf_file=f'{cif_path_without_ext}_opt_n.gjf', osc_type='n')
    else:
        raise OSCTypeError

    skip_specified_cal = []
    if args.read:
        print('Skip calculation of reorganization energy.')
    elif args.resume:
        rcal.check_extension_log(f'{cif_path_without_ext}_opt_n.gjf')
        skip_specified_cal = check_reorganization_energy_completion(cif_path_without_ext, args.osc_type, extension_log=rcal._extension_log)
    else:
        print('Calculate reorganization energy.')

    reorg_energy = rcal.calc_reorganization(gau_com=gau_com, only_read=args.read, is_output_detail=True, skip_specified_cal=skip_specified_cal)

    print_reorg_energy(args.osc_type, reorg_energy)

    ##### Calculate transfer integrals #####
    transfer_integrals = []
    mom_dis_ti = [] # Store moment of inertia, distance between centers of weight and transfer integral

    expand_mols = cif_reader.expand_mols(args.cellsize)
    for s in range(len(cif_reader.unique_symbols.keys())):
        unique_symbols = cif_reader.unique_symbols[s]
        unique_coords = cif_reader.unique_coords[s]
        unique_coords = cif_reader.convert_frac_to_cart(unique_coords)
        for (i, j, k), expand_mol in expand_mols.items():
            for t, (symbols, coordinates) in expand_mol.items():
                # Skip creating gjf for transfer integrals because they are molecules with translation symmetry
                if s > t:
                    continue
                elif s == t:
                    if (i, j, k) == (0, 0, 0):
                        continue
                    elif i < 0 or (i == 0 and (j < 0 or (j == 0 and k < 0))):
                        continue

                coordinates = cif_reader.convert_frac_to_cart(coordinates)

                min_distance = cal_min_distance(
                    unique_symbols, unique_coords,
                    symbols, coordinates
                )
                if min_distance > 5:
                    print()
                    print(f'Skip calculation of transfer integral from {s}-th in (0,0,0) cell to {t}-th in ({i},{j},{k}) cell because the minimum distance is over 5 \u212B.\n')
                    continue

                moment, _ = cal_moment_of_inertia(
                    unique_symbols, unique_coords,
                    symbols, coordinates
                )

                distance = cal_distance_between_cen_of_weight(
                    unique_symbols, unique_coords,
                    symbols, coordinates
                )

                is_run_ti = True
                same_ti = 0

                # skip calculation of transfer integrals using moment of inertia and distance between centers of weight.
                if not args.fullcal:
                    for m, d, ti in mom_dis_ti:
                        if (np.all(m - MOMENT_OF_INERTIA_ERROR < moment) and np.all(moment < m + MOMENT_OF_INERTIA_ERROR)) and (d - CENTER_OF_WEIGHT_ERROR < distance < d + CENTER_OF_WEIGHT_ERROR):
                            is_run_ti = False
                            same_ti = ti
                            break

                if is_run_ti:
                    gjf_name = f'{filename}-({s}_{t}_{i}_{j}_{k})'
                    gjf_file = f'{directory}/{gjf_name}'

                    tcal = Tcal(gjf_file)

                    is_normal_term = False
                    if args.resume:
                        tcal.check_extension_log()
                        is_normal_term = check_transfer_integral_completion(gjf_file, extension_log=tcal._extension_log)

                    if not args.read and not is_normal_term:
                        print()
                        print('Create gjf for transfer integral.')
                        create_ti_gjf(
                            {'symbols': unique_symbols, 'coordinates': unique_coords},
                            {'symbols': symbols, 'coordinates': coordinates},
                            gjf_basename=gjf_name,
                            save_dir=directory,
                            cpu=args.cpu,
                            mem=args.mem,
                            method=args.method,
                        )
                        tcal.create_monomer_file()

                        if args.g09:
                            gaussian_command = 'g09'
                        else:
                            gaussian_command = 'g16'
                        print(f'Calculate transfer integral from {s}-th in (0,0,0) cell to {t}-th in ({i},{j},{k}) cell.')
                        tcal.run_gaussian(gaussian_command)
                    else:
                        print()
                        print(f'Skip calculation of transfer integral from {s}-th in (0,0,0) cell to {t}-th in ({i},{j},{k}) cell.')

                    tcal.check_extension_log()
                    tcal.read_monomer1()
                    tcal.read_monomer2()
                    tcal.read_dimer()

                    if args.osc_type == 'p':
                        transfer = Tcal.cal_transfer_integrals(
                            tcal.mo1[tcal.n_elect1-1], tcal.overlap, tcal.fock, tcal.mo2[tcal.n_elect2-1]
                        )
                    elif args.osc_type == 'n':
                        transfer = Tcal.cal_transfer_integrals(
                            tcal.mo1[tcal.n_elect1], tcal.overlap, tcal.fock, tcal.mo2[tcal.n_elect2]
                        )

                    transfer = transfer * 1e-3 # meV to eV
                    print_transfer_integral(args.osc_type, transfer)
                    transfer_integrals.append((s, t, i, j, k, transfer))
                    mom_dis_ti.append((moment, distance, transfer))
                else:
                    print()
                    print(f'Skip calculation of transfer integral from {s}-th in (0,0,0) cell to {t}-th in ({i},{j},{k}) cell due to identical moment of inertia and distance between centers of weight.')
                    print_transfer_integral(args.osc_type, same_ti)
                    transfer_integrals.append((s, t, i, j, k, same_ti))

    ##### Calculate mobility tensor considering anisotropy. #####
    hop = []

    for s, t, i, j, k, ti in transfer_integrals:
        hop.append((s, t, i, j, k, marcus_rate(ti, reorg_energy)))

    diffusion_coef_tensor = diffusion_coefficient_tensor(cif_reader.lattice * 1e-8, hop)
    print_tensor(diffusion_coef_tensor, msg="Diffusion coefficient tensor (cm^2/s)")
    mu = mobility_tensor(diffusion_coef_tensor)
    print_tensor(mu, msg="Mobility tensor (cm^2/Vs)")
    value, vector = cal_eigenvalue_decomposition(mu)
    print_mobility(value, vector)

    ##### Simulate mobility tensor calculation using Monte Carlo method #####
    if args.mc:
        D_MC = diffusion_coefficient_tensor_MC(cif_reader.lattice * 1e-8, hop)
        print_tensor(D_MC, msg="Diffusion coefficient tensor (cm^2/s) (MC)")
        mu_MC = mobility_tensor(D_MC)
        print_tensor(mu_MC, msg="Mobility tensor (cm^2/Vs) (MC)")
        value_MC, vector_MC = cal_eigenvalue_decomposition(mu_MC)
        print_mobility(value_MC, vector_MC, sim_type='MC')

    ##### Simulate mobility tensor calculation using Ordinary Differential Equation method #####
    if args.ode:
        D_ODE = diffusion_coefficient_tensor_ODE(cif_reader.lattice * 1e-8, hop)
        print_tensor(D_ODE, msg="Diffusion coefficient tensor (cm^2/s) (ODE)")
        mu_ODE = mobility_tensor(D_ODE)
        print_tensor(mu_ODE, msg="Mobility tensor (cm^2/Vs) (ODE)")
        value_ODE, vector_ODE = cal_eigenvalue_decomposition(mu_ODE)
        print_mobility(value_ODE, vector_ODE, sim_type='ODE')

    # Save reorganization, transfer integrals, hop, mobility tensor
    if args.pickle:
        with open(f'{cif_path_without_ext}_result.pkl', 'wb') as f:
            pickle.dump({
                'osc_type': args.osc_type,
                'lattice': cif_reader.lattice,
                'z_value': cif_reader.z_value,
                'reorganization': reorg_energy,
                'transfer_integrals': transfer_integrals,
                'hop': hop,
                'diffusion_coefficient_tensor': diffusion_coef_tensor,
                'mobility_tensor': mu,
                'mobility_value': value,
                'mobility_vector': vector
            }, f)

    Tcal.print_timestamp()
    end_time = time()
    elapsed_time = end_time - start_time
    elapsed_time_h = int(elapsed_time // 3600)
    elapsed_time_min = int((elapsed_time - elapsed_time_h * 3600) // 60)
    elapsed_time_sec = int(elapsed_time - elapsed_time_h * 3600 - elapsed_time_min * 60)
    elapsed_time_ms = (elapsed_time - elapsed_time_h * 3600 - elapsed_time_min * 60 - elapsed_time_sec) * 1000
    if elapsed_time < 1:
        print(f'Elapsed Time: {elapsed_time_ms:.0f} ms')
    elif elapsed_time < 60:
        print(f'Elapsed Time: {elapsed_time_sec} sec')
    elif elapsed_time < 3600:
        print(f'Elapsed Time: {elapsed_time_min} min {elapsed_time_sec} sec')
    else:
        print(f'Elapsed Time: {elapsed_time_h} h {elapsed_time_min} min {elapsed_time_sec} sec')


def atom_weight(symbol: str) -> float:
    """Get atom weight

    Parameters
    ----------
    symbol : str
        Symbol of atom

    Returns
    -------
    float
        Atomic weight
    """
    ELEMENT_PROP = CifReader.ELEMENT_PROP
    weight = ELEMENT_PROP[ELEMENT_PROP['symbol'] == symbol]['weight'].values[0]

    return weight


def cal_cen_of_weight(
    symbols1: NDArray[str],
    coordinates1: NDArray[np.float64],
    symbols2: Optional[NDArray[str]] = None,
    coordinates2: Optional[NDArray[np.float64]] = None,
) -> NDArray[np.float64]:
    """Calculate center of weight

    Parameters
    ----------
    symbols1 : NDArray[str]
        Symbols of atoms in one monomer
    coordinates1 : NDArray[np.float64]
        Coordinates of atoms in one monomer
    symbols2 : Optional[NDArray[str]], optional
        Symbols of atoms in another monomer, by default None
    coordinates2 : Optional[NDArray[np.float64]], optional
        Coordinates of atoms in another monomer, by default None

    Returns
    -------
    NDArray[np.float64]
        Center of weight
    """
    if symbols2 is not None and coordinates2 is not None:
        symbols1 = np.concatenate((symbols1, symbols2), axis=0)
        coordinates1 = np.concatenate((coordinates1, coordinates2), axis=0)

    weights = np.array([atom_weight(sym) for sym in symbols1])
    total_weight = np.sum(weights)

    weighted_coords = weights[:, np.newaxis] * coordinates1
    weighted_sum = np.sum(weighted_coords, axis=0)

    cen_of_weight = weighted_sum / total_weight

    return cen_of_weight


def cal_distance_between_cen_of_weight(
    symbols1: NDArray[str],
    coordinates1: NDArray[np.float64],
    symbols2: NDArray[str],
    coordinates2: NDArray[np.float64],
) -> float:
    """Calculate distance between centers of weight

    Parameters
    ----------
    symbols1 : NDArray[str]
        Symbols of atoms in one monomer
    coordinates1 : NDArray[np.float64]
        Coordinates of atoms in one monomer
    symbols2 : NDArray[str]
        Symbols of atoms in another monomer
    coordinates2 : NDArray[np.float64]
        Coordinates of atoms in another monomer

    Returns
    -------
    float
        Distance between centers of weight
    """
    mol1_cen_coord = cal_cen_of_weight(symbols1, coordinates1)
    mol2_cen_coord = cal_cen_of_weight(symbols2, coordinates2)
    distance = np.sqrt(np.sum(np.square(mol1_cen_coord-mol2_cen_coord)))

    return distance


def cal_eigenvalue_decomposition(mobility_tensor: NDArray[np.float64]) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Calculate eigenvalue decomposition of mobility tensor

    Parameters
    ----------
    mobility_tensor : NDArray[np.float64]
        Mobility tensor

    Returns
    -------
    Tuple[NDArray[np.float64], NDArray[np.float64]]
        Eigenvalue(mobility value) and eigenvector(mobility vector)
    """
    value, vector = np.linalg.eig(mobility_tensor)
    return value, vector


def cal_min_distance(
    symbols1: NDArray[str],
    coords1: NDArray[np.float64],
    symbols2: NDArray[str],
    coords2: NDArray[np.float64],
) -> float:
    """Calculate minimum distance between two sets of atoms.

    Parameters
    ----------
    symbols1 : NDArray[str]
        Symbols of atoms in one monomer
    coords1 : NDArray[np.float64]
        Coordinates of atoms in one monomer
    symbols2 : NDArray[str]
        Symbols of atoms in another monomer
    coords2 : NDArray[np.float64]
        Coordinates of atoms in another monomer

    Returns
    -------
    float
        Minimum distance between two sets of atoms
    """
    ELEMENT_PROP = CifReader.ELEMENT_PROP
    VDW_RADII = ELEMENT_PROP[['symbol', 'vdw_radius']].set_index('symbol').to_dict()['vdw_radius']

    radii1 = np.array(
        [VDW_RADII[symbol] for symbol in symbols1]
    )
    radii2 = np.array(
        [VDW_RADII[symbol] for symbol in symbols2]
    )

    distances = np.sqrt(np.sum((coords1[:, np.newaxis] - coords2)**2, axis=2)) - radii1[:, np.newaxis] - radii2

    min_distance = np.min(distances)

    return min_distance


def cal_moment_of_inertia(
    symbols1: NDArray[str],
    coordinates1: NDArray[np.float64],
    symbols2: NDArray[str],
    coordinates2: NDArray[np.float64],
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Calculate moment of inertia and eigenvectors of the inertia tensor.

    Parameters
    ----------
    symbols1 : NDArray[str]
        Symbols of atoms in one monomer
    coordinates1 : NDArray[np.float64]
        Coordinates of atoms in one monomer
    symbols2 : NDArray[str]
        Symbols of atoms in another monomer
    coordinates2 : NDArray[np.float64]
        Coordinates of atoms in another monomer

    Returns
    -------
    Tuple[NDArray[np.float64], NDArray[np.float64]]
        Moment of inertia and eigenvectors of the inertia tensor
    """
    symbols1 = np.concatenate((symbols1, symbols2), axis=0)
    coordinates1 = np.concatenate((coordinates1, coordinates2), axis=0)

    cen_of_weight = cal_cen_of_weight(symbols1, coordinates1)

    weights = np.array([atom_weight(sym) for sym in symbols1])

    xi = coordinates1[:, 0] - cen_of_weight[0]
    yi = coordinates1[:, 1] - cen_of_weight[1]
    zi = coordinates1[:, 2] - cen_of_weight[2]

    tmp_coords = np.column_stack((xi, yi, zi))

    moment = np.zeros((3, 3))

    for i in range(3):
        moment[i, i] = np.sum(weights * (tmp_coords[:, (i+1)%3]**2 + tmp_coords[:, (i+2)%3]**2))

    for i in range(3):
        for j in range(i+1, 3):
            moment[i, j] = moment[j, i] = -np.sum(weights * tmp_coords[:, i] * tmp_coords[:, j])

    moment, p = np.linalg.eig(moment)

    return moment, p


def check_reorganization_energy_completion(
    cif_path_without_ext: str,
    osc_type: Literal['p', 'n'],
    extension_log: str = '.log'
) -> List[Literal['opt_neutral', 'opt_ion', 'neutral', 'ion']]:
    """Check if all reorganization energy calculations are completed normally.

    Parameters
    ----------
    cif_path_without_ext : str
        Base path of cif file (without extension)
    osc_type : Literal['p', 'n']
        Semiconductor type (p-type or n-type)
    extension_log : str
        Extension of log file

    Returns
    -------
    List[Literal['opt_neutral', 'opt_ion', 'neutral', 'ion']]
        List of calculations to skip
    """
    skip_specified_cal = []
    if check_normal_termination(f'{cif_path_without_ext}_opt_n{extension_log}'):
        skip_specified_cal.append('opt_neutral')
    if check_normal_termination(f'{cif_path_without_ext}_n{extension_log}'):
        skip_specified_cal.append('neutral')

    if osc_type == 'p':
        if check_normal_termination(f'{cif_path_without_ext}_opt_c{extension_log}'):
            skip_specified_cal.append('opt_ion')
        if check_normal_termination(f'{cif_path_without_ext}_c{extension_log}'):
            skip_specified_cal.append('ion')
    elif osc_type == 'n':
        if check_normal_termination(f'{cif_path_without_ext}_opt_a{extension_log}'):
            skip_specified_cal.append('opt_ion')
        if check_normal_termination(f'{cif_path_without_ext}_a{extension_log}'):
            skip_specified_cal.append('ion')

    return skip_specified_cal


def check_transfer_integral_completion(gjf_file: str, extension_log: str = '.log') -> bool:
    """Check if all transfer integral calculations are completed normally.

    Parameters
    ----------
    gjf_file : str
        Base path of gjf file (without extension)

    Returns
    -------
    bool
        True if all calculations (dimer, monomer1, monomer2) terminated normally
    """
    required_files = ['', '_m1', '_m2']
    return all(
        check_normal_termination(f'{gjf_file}{suffix}{extension_log}')
        for suffix in required_files
    )


def create_reorg_gjf(
    symbols: NDArray[str],
    coordinates: NDArray[np.float64],
    basename: str,
    save_dir: str,
    cpu: int,
    mem: int,
    method: str,
) -> None:
    """Create gjf file for reorganization energy calculation.

    Parameters
    ----------
    symbols : NDArray[str]
        Symbols of atoms
    coordinates : NDArray[np.float64]
        Coordinates of atoms
    basename : str
        Base name of gjf file
    save_dir : str
        Directory to save gjf file
    cpu : int
        Number of cpu
    mem : int
        Number of memory [GB]
    method : str
        Calculation method used in Gaussian calculations
    """
    gjf_maker = GjfMaker()
    gjf_maker.set_function(method)
    gjf_maker.create_chk_file()
    gjf_maker.output_detail()
    gjf_maker.opt()

    gjf_maker.set_symbols(symbols)
    gjf_maker.set_coordinates(coordinates)
    gjf_maker.set_resource(cpu_num=cpu, mem_num=mem)

    gjf_maker.export_gjf(
        file_name=f'{basename}_opt_n',
        save_dir=save_dir,
        chk_rwf_name=f'{save_dir}/{basename}_opt_n'
    )


def create_ti_gjf(
    unique_mol: Dict[str, Union[NDArray[str], NDArray[np.float64]]],
    neighbor_mol: Dict[str, Union[NDArray[str], NDArray[np.float64]]],
    gjf_basename: str,
    save_dir: str = '.',
    cpu: int = 4,
    mem: int = 16,
    method: str = 'B3LYP/6-31G*',
) -> None:
    """Create gjf file for transfer integral calculation.

    Parameters
    ----------
    unique_mol : Dict[str, Union[NDArray[str], NDArray[np.float64]]]
        Dictionary containing symbols and coordinates of unique monomer
    neighbor_mol : Dict[str, Union[NDArray[str], NDArray[np.float64]]]
        Dictionary containing symbols and coordinates of neighbor monomer
    gjf_basename : str
        Base name of gjf file
    save_dir : str
        Directory to save gjf file, by default '.'
    cpu : int
        Number of cpu, by default 4
    mem : int
        Number of memory [GB], by default 16
    method : str
        Calculation method used in Gaussian calculations, by default 'B3LYP/6-31G(d,p)'
    """
    gjf_maker = GjfMaker()
    gjf_maker.set_resource(cpu_num=cpu, mem_num=mem)
    gjf_maker.set_function(method)
    gjf_maker.create_chk_file()
    gjf_maker.add_root('Symmetry=None')

    gjf_maker.set_symbols(unique_mol['symbols'])
    gjf_maker.set_coordinates(unique_mol['coordinates'])
    gjf_maker.set_symbols(neighbor_mol['symbols'])
    gjf_maker.set_coordinates(neighbor_mol['coordinates'])

    gjf_maker.add_link()
    gjf_maker.add_root('Symmetry=None')
    gjf_maker.add_root('Pop=Full')
    gjf_maker.add_root('IOp(3/33=4,5/33=3)')

    gjf_maker.export_gjf(file_name=gjf_basename, save_dir=save_dir)


def print_mobility(value: NDArray[np.float64], vector: NDArray[np.float64], sim_type: Literal['MC', 'ODE'] = ''):
    """Print mobility and mobility vector

    Parameters
    ----------
    value : NDArray[np.float64]
        Mobility value
    vector : NDArray[np.float64]
        Mobility vector
    sim_type : str
        Simulation type (MC or ODE)
    """
    msg_value = 'Mobility eigenvalues (cm^2/Vs)'
    msg_vector = 'Mobility eigenvectors'
    direction = ['x', 'y', 'z']

    if sim_type:
        msg_value += f' ({sim_type})'
        msg_vector += f' ({sim_type})'

    print()
    print('-' * (len(msg_value)+2))
    print(f' {msg_value} ')
    print('-' * (len(msg_value)+2))
    print(f"{value[0]:12.6g} {value[1]:12.6g} {value[2]:12.6g}")
    print()

    print()
    print('-' * (len(msg_vector)+2))
    print(f' {msg_vector} ')
    print('-' * (len(msg_vector)+2))
    print('       vector1      vector2      vector3')
    for v, d in zip(vector, direction):
        print(f'{d} {v[0]:12.6g} {v[1]:12.6g} {v[2]:12.6g}')
    print()


def print_reorg_energy(osc_type: Literal['p', 'n'], reorg_energy: float):
    """Print reorganization energy

    Parameters
    ----------
    osc_type : Literal['p', 'n']
        Semiconductor type (p-type or n-type)
    reorg_energy : float
        Reorganization energy [eV]
    """
    print()
    print('-----------------------')
    print(' Reorganization energy ')
    print('-----------------------')
    print(f'{osc_type}-type: {reorg_energy:10.6g} eV\n')


def print_tensor(mu: NDArray[np.float64], msg: str = 'Mobility tensor'):
    """Print mobility tensor

    Parameters
    ----------
    mu : NDArray[np.float64]
        Mobility tensor
    msg : str
        Message, by default 'Mobility tensor'
    """
    print()
    print('-' * (len(msg)+2))
    print(f' {msg} ')
    print('-' * (len(msg)+2))
    for a in mu:
        print(f"{a[0]:12.6g} {a[1]:12.6g} {a[2]:12.6g}")
    print()


def print_transfer_integral(osc_type: Literal['p', 'n'], transfer: float):
    """Print transfer integral

    Parameters
    ----------
    osc_type : Literal['p', 'n']
        Semiconductor type (p-type or n-type)
    transfer : float
        Transfer integral [eV]
    """
    mol_orb = {'p': 'HOMO', 'n': 'LUMO'}
    print()
    print('-------------------')
    print(' Transfer integral ')
    print('-------------------')
    print(f'{mol_orb[osc_type]}: {transfer:12.6g} eV\n')


def read_pickle(file_name: str):
    print(f'\nInput File Name: {file_name}')

    with open(file_name, 'rb') as f:
        results = pickle.load(f)

    # print(results)

    print(f'\nCalculate as {results["osc_type"]}-type organic semiconductor.')

    print_reorg_energy(results['osc_type'], results['reorganization'])

    for s, t, i, j, k, ti in results['transfer_integrals']:
        print()
        print(f'{s}-th in (0,0,0) cell to {t}-th in ({i},{j},{k}) cell')
        print_transfer_integral(results['osc_type'], ti)

    print_tensor(results['diffusion_coefficient_tensor'], msg="Diffusion coefficient tensor (cm^2/s)")

    print_tensor(results['mobility_tensor'], msg="Mobility tensor (cm^2/Vs)")

    print_mobility(results['mobility_value'], results['mobility_vector'])


class OSCTypeError(Exception):
    """Exception for semiconductor type"""
    pass


if __name__ == '__main__':
    main()
