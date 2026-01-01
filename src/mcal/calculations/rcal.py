"""Rcal"""
import argparse
import functools
import os
import subprocess
from datetime import datetime
from pathlib import Path
from time import time
from typing import List, Literal

from mcal.utils.cif_reader import CifReader
from mcal.utils.gjf_maker import GjfMaker


print = functools.partial(print, flush=True)


def main():
    """This code is to execute rcal for command line."""
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='cif file name or gjf file name', type=str)
    parser.add_argument('osc_type', help='organic semiconductor type', type=str)
    parser.add_argument(
        '-M', '--method',
        help='calculation method used in Gaussian calculations (default is B3LYP/6-31G(d,p)). ' \
            + 'But if you use a gjf file instead of a cif file, the method in the gjf file will be used',
        type=str,
        default='B3LYP/6-31G(d,p)',
    )
    parser.add_argument(
        '-c', '--cpu',
        help='setting the number of cpu (default is 4). ' \
            + 'But if you use a gjf file instead of a cif file, the number of cpu in the gjf file will be used',
        type=int, default=4
    )
    parser.add_argument(
        '-m', '--mem',
        help='setting the number of memory [GB] (default is 10 GB). ' \
            + 'But if you use a gjf file instead of a cif file, the number of memory in the gjf file will be used',
        type=int,
        default=10,
    )
    parser.add_argument('-g', '--g09', help='use Gaussian 09 (default is Gaussian 16)', action='store_true')
    parser.add_argument('-r', '--read', help='read log files without executing Gaussian', action='store_true')
    args = parser.parse_args()

    print('---------------------------------------')
    print(' rcal beta (2025/06/21) by Matsui Lab. ')
    print('---------------------------------------')
    print(f'\nInput File Name: {args.file}')
    Rcal.print_timestamp()
    before = time()

    if args.file.endswith('.cif'):
        cif_file = Path(args.file)
        directory = cif_file.parent
        filename = cif_file.stem.replace('_opt_n', '')
        file_path_without_ext = f'{directory}/{filename}_opt_n'

        cif_reader = CifReader(cif_path=cif_file)
        symbols = cif_reader.unique_symbols[0]
        coordinates = cif_reader.unique_coords[0]
        coordinates = cif_reader.convert_frac_to_cart(coordinates)

        gjf_maker = GjfMaker()
        gjf_maker.create_chk_file()
        gjf_maker.output_detail()
        gjf_maker.opt()

        gjf_maker.set_symbols(symbols)
        gjf_maker.set_coordinates(coordinates)
        gjf_maker.set_function(args.method)
        gjf_maker.set_charge_spin(charge=0, spin=1)
        gjf_maker.set_resource(cpu_num=args.cpu, mem_num=args.mem)

        gjf_maker.export_gjf(
            file_name=f'{file_path_without_ext}',
            chk_rwf_name=f'{file_path_without_ext}',
        )
    elif args.file.endswith('.gjf'):
        gjf_file = Path(args.file)
        directory = gjf_file.parent
        filename = gjf_file.stem

        file_path_without_ext = f'{directory}/{filename}'
    else:
        raise ValueError('Input file must be a cif file or a gjf file.')

    if args.osc_type.lower() == 'p':
        rcal = Rcal(gjf_file=f'{file_path_without_ext}.gjf')
    elif args.osc_type.lower() == 'n':
        rcal = Rcal(gjf_file=f'{file_path_without_ext}.gjf', osc_type='n')
    else:
        raise OSCTypeError

    if args.g09:
        gau_com = 'g09'
    else:
        gau_com = 'g16'

    reorg_energy = rcal.calc_reorganization(gau_com=gau_com, only_read=args.read, is_output_detail=True)

    print()
    print('-----------------------')
    print(' Reorganization energy ')
    print('-----------------------')
    print(f'{reorg_energy:12.6f} eV')
    print()

    Rcal.print_timestamp()
    after = time()
    print(f'Elapsed Time: {(after - before) * 1000:.0f} ms')


class Rcal:
    """Calculate organization energy."""
    def __init__(self, gjf_file: str, osc_type: Literal['p', 'n'] = 'p'):
        """
        Initialize Rcal.

        Parameters
        ----------
        gjf_file : str
            gjf file name.
        osc_type : Literal['p', 'n']
            organic semiconductor type, 'p' is positive, 'n' is negative, by default 'p'.
        """
        self.gjf_file = gjf_file
        self.ion = None
        self._extension_log = '.log'
        self._gjf_lines = {'%': [], '#': []}

        self._input_gjf()

        if osc_type.lower() == 'p':
            self.ion = 'c'
        elif osc_type.lower() == 'n':
            self.ion = 'a'

    @ staticmethod
    def check_error_term(line: str) -> None:
        """
        Check the error term of Gaussian.

        Parameters
        ----------
        line : str
            last line of the log file.

        Raises
        ------
        GausTermError
            if the calculation of Gaussian was error termination.
        """
        line = line.strip()

        if 'Normal termination' not in line:
            raise GausTermError('The calculation of Gaussian was error termination.')

    @staticmethod
    def print_timestamp() -> None:
        """Print timestamp."""
        month = {
            1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
            7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec',
        }
        dt_now = datetime.now()
        print(f"Timestamp: {dt_now.strftime('%a')} {month[dt_now.month]} {dt_now.strftime('%d %H:%M:%S %Y')}")

    def calc_reorganization(
        self,
        gau_com: str = 'g16',
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_specified_cal: List[Literal['opt_neutral', 'opt_ion', 'neutral', 'ion']] = [],
    ) -> float:
        """
        Calculate reorganization energy.

        Parameters
        ----------
        gau_com : str
            Gaussian command, by default 'g16'.
        only_read : bool
            if True, the calculation is only read, by default False.
        is_output_detail : bool
            if True, the calculation detail will be output, by default False.
        skip_specified_cal : List[Literal['opt_neutral', 'opt_ion', 'neutral', 'ion']]
            if specified, the calculation of the specified type will be skipped, by default [].

        Returns
        -------
        float
            reorganization energy [eV].
        """
        file_path = Path(self.gjf_file)
        filename = file_path.stem.replace('_opt_n', '')
        directory = file_path.parent
        basename = f'{directory}/{filename}'

        energy = []

        # 中性分子の構造最適化とエネルギー計算
        only_read_opt_n = only_read
        if not only_read and 'opt_neutral' not in skip_specified_cal:
            print('>', gau_com, self.gjf_file)
            subprocess.run([gau_com, self.gjf_file])

        skip_opt_neutral = True if 'opt_neutral' in skip_specified_cal else False

        energy.append(self.extract_energy(self.gjf_file, only_read=only_read_opt_n, is_output_detail=is_output_detail, skip_cal=skip_opt_neutral))

        # カチオンかアニオンのエネルギー計算
        only_read_ion = only_read
        previous_name, _ = os.path.splitext(self.gjf_file)
        gjf = f'{basename}_{self.ion}.gjf'
        if not only_read and 'ion' not in skip_specified_cal:
            self._create_gjf(file_name=gjf, prevous_name=previous_name, ion=self.ion)
            print('>', gau_com, gjf)
            subprocess.run([gau_com, gjf])

        skip_ion = True if 'ion' in skip_specified_cal else False

        energy.append(self.extract_energy(gjf, only_read=only_read_ion, is_output_detail=is_output_detail, skip_cal=skip_ion))

        # カチオンかアニオンの構造最適化とエネルギー計算
        only_read_opt_ion = only_read
        gjf = f'{basename}_opt_{self.ion}.gjf'
        if not only_read and 'opt_ion' not in skip_specified_cal:
            self._create_gjf(file_name=gjf, prevous_name=previous_name, ion=self.ion, is_opt=True)
            print('>', gau_com, gjf)
            subprocess.run([gau_com, gjf])

        skip_opt_ion = True if 'opt_ion' in skip_specified_cal else False

        energy.append(self.extract_energy(gjf, only_read=only_read_opt_ion, is_output_detail=is_output_detail, skip_cal=skip_opt_ion))


        # 中性分子のエネルギー計算
        only_read_neutral = only_read
        previous_name, _ = os.path.splitext(gjf)
        ion = 'n'
        gjf = f'{basename}_{ion}.gjf'
        if not only_read and 'neutral' not in skip_specified_cal:
            self._create_gjf(file_name=gjf, prevous_name=previous_name, ion=ion)
            print('>', gau_com, gjf)
            subprocess.run([gau_com, gjf])

        skip_neutral = True if 'neutral' in skip_specified_cal else False

        energy.append(self.extract_energy(gjf, only_read=only_read_neutral, is_output_detail=is_output_detail, skip_cal=skip_neutral))

        return ((energy[3] - energy[2]) + (energy[1] - energy[0]))

    def check_extension_log(self, gjf: str) -> None:
        """Check the extension of log file.

        Parameters
        ----------
        gjf : str
            gjf file name.
        """
        if os.path.exists(f'{os.path.splitext(gjf)[0]}.out'):
            self._extension_log = '.out'
        else :
            self._extension_log = '.log'

    def extract_energy(
        self,
        gjf: str,
        only_read: bool = False,
        is_output_detail: bool = False,
        skip_cal: bool = False
    ) -> float:
        """Extract energy from log file.

        Parameters
        ----------
        gjf : str
            gjf file name.
        only_read : bool
            if True, the calculation is only read, by default False.
        is_output_detail : bool
            if True, the calculation detail will be output, by default False.

        Returns
        -------
        float
            total energy.
        """
        self.check_extension_log(gjf)
        log_file = f'{os.path.splitext(gjf)[0]}{self._extension_log}'

        with open(log_file) as f:
            last_line = ''
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip()

                if line:
                    last_line = line

                if line.startswith('SCF Done:'):
                    energy = float(line.split()[4]) * 27.2114

            self.check_error_term(last_line)

            if is_output_detail:
                gjf = Path(gjf)
                if not only_read and not skip_cal:
                    print(f'{gjf} calculation completed.')
                elif skip_cal:
                    print(f'{gjf} calculation skipped.')

                print(f'reading {gjf.parent}/{gjf.stem}{self._extension_log}')
                print()
                print('--------------')
                print(' Total energy ')
                print('--------------')
                print(f'{energy:12.6f} eV')
                print()

            return energy

    def _create_gjf(
        self,
        file_name: str,
        prevous_name: str,
        ion: Literal['c', 'a', 'n'],
        is_opt: bool = False,
    ) -> None:
        """
        Create gjf file.

        Parameters
        ----------
        file_name : str
            file name.
        prevous_name : str
            previous file name.
        ion : Literal['c', 'a', 'n']
            ion type. 'c' is cation, 'a' is anion, 'n' is neutral molecule.
        is_opt : bool
            if True, the calculation is optimization, by default False.
        """
        file_name, _ = os.path.splitext(file_name)

        with open(f'{file_name}.gjf', 'w') as f:
            for line in self._gjf_lines['%']:
                if r'%oldchk' in line.lower():
                    continue
                elif r'%chk' in line.lower():
                    continue
                else:
                    f.write(line)

            f.write(f'%oldchk={prevous_name}.chk\n')
            if is_opt:
                f.write(f'%chk={file_name}.chk\n')

            for line in self._gjf_lines['#']:
                if 'geom' in line.lower():
                    continue
                elif 'opt' in line.lower():
                    continue
                else:
                    f.write(line)

            f.write('#  Geom=Checkpoint\n')
            if is_opt:
                f.write('#  Opt=Tight\n')
            f.write('\n')
            f.write('Defalut Title\n')
            f.write('\n')

            if ion == 'c':
                f.write('1 2\n\n')
            elif ion == 'a':
                f.write('-1, 2\n\n')
            else:
                f.write('0 1\n\n')

    def _input_gjf(self) -> None:
        """Input link 0 command and root options from gjf file."""
        with open(self.gjf_file, 'r') as f:
            for line in f:
                if line.startswith('%'):
                    self._gjf_lines['%'].append(line)
                elif line.startswith('#'):
                    self._gjf_lines['#'].append(line)
                elif 'link' in line.lower():
                    raise ValueError("Please do not use Link.")


class GausTermError(Exception):
    """Exception for Gaussian error termination."""
    pass


class OSCTypeError(Exception):
    """Exception for organic semiconductor type."""
    pass


if __name__ == '__main__':
    main()
