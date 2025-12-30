""" GjfMaker beta (2025/08/18)"""
import os
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


class GjfMaker:
    parent_dir = Path(os.path.abspath(__file__)).parent.parent
    ELEMENT_PROP = pd.read_csv(f'{parent_dir}/constants/element_properties.csv')
    ELEMENTS_NUM = ELEMENT_PROP[['symbol', 'number']].set_index('symbol').to_dict()['number']

    def __init__(self, remove_radical_flag: bool = False) -> None:
        """Initialize the GjfMaker class.

        Parameters
        ----------
        remove_radical_flag : bool
            If True, the radical is removed, by default False.
            Not recommended to use this flag.
        """
        self.remove_radical_flag = remove_radical_flag

        self._link_dict = {}
        self._link_num = 0
        self._cpu_num = '%NProcShared=9\n'
        self._mem_num = '%Mem=40GB\n'
        self._is_chk_file = False
        self._is_rwf_file = False
        self.is_detail = False
        self._function = '# B3LYP/6-31G(d,p)\n'
        self._root = []
        self._title = '\nDefalut Title\n'
        self._symbols = []
        self._charge_spin = ''
        self._coordinates = []

        self._root_space = ' '

    def add_link(self) -> None:
        """Add the link of the gjf file."""
        self.create_chk_file()
        self._link_dict[self._link_num] = [self._cpu_num, self._mem_num, self._function, self._root]
        self._link_num += 1
        self._root = ['Geom=AllCheck\n', 'Guess=Read\n']

    def add_root(self, root: str) -> None:
        """Add the root of the gjf file."""
        if root[0] == '#':
            root = root[1:]
        self._root.append(root + '\n')

    def create_chk_file(self) -> None:
        """Create the chk file."""
        self._is_chk_file = True

    def create_rwf_file(self) -> None:
        """Create the rwf file."""
        self._is_rwf_file = True

    def export_gjf(
        self,
        file_name: str,
        save_dir: str = '.',
        chk_rwf_name: Optional[str] = None
    ) -> None:
        """Export the gjf file.

        Parameters
        ----------
        file_name : str
            Name of the gjf file.
        save_dir : str
            Directory to save the gjf file.
        chk_rwf_name : Optional[str]
            Name of the chk or rwf file.
        """
        if self.remove_radical_flag and self._check_radical():
            print('This molecule is radical')
            return

        os.makedirs(save_dir, exist_ok=True)

        self._link_dict[self._link_num] = [self._cpu_num, self._mem_num, self._function, self._root]

        if self._is_chk_file:
            if chk_rwf_name:
                self._chk_file = f'%Chk={chk_rwf_name}.chk\n'
            else:
                self._chk_file = f'%Chk={save_dir}/{file_name}.chk\n'
        else:
            self._chk_file = ''

        if self._is_rwf_file:
            if chk_rwf_name:
                self._rwf_file = f'%RWF={chk_rwf_name}.rwf\n'
            else:
                self._rwf_file = f'%RWF={save_dir}/{file_name}.rwf\n'
        else:
            self._rwf_file = ''

        with open(f'{save_dir}/{file_name}.gjf', 'w') as f:
            for i in range(self._link_num+1):
                _cpn_num, _mem_num, _function, _root = self._link_dict[i]
                if i > 0:
                    f.write(f'\n--Link{i}--\n')
                f.write(_cpn_num)
                f.write(_mem_num)
                f.write(self._chk_file)
                f.write(self._rwf_file)
                f.write(_function)
                for r in _root:
                    f.write(f'#{self._root_space}{r}')

                if i == 0:
                    f.write(self._title)

                    if not self._charge_spin:
                        self._check_charge_spin()
                    f.write(self._charge_spin)

                    for sym, (x, y, z) in zip(self._symbols, self._coordinates):
                        x = f'{x:.8f}'.rjust(14, ' ')
                        y = f'{y:.8f}'.rjust(14, ' ')
                        z = f'{z:.8f}'.rjust(14, ' ')
                        f.write(f' {sym}     {x} {y} {z}\n')
            f.write('\n')

    def opt(self, tight: bool = True) -> None:
        """Set the optimization option of the gjf file.

        Parameters
        ----------
        tight : bool
            If True, the optimization is tight, by default True.
        """
        opt = 'Opt'
        if tight:
            opt = 'Opt=Tight'
        self.add_root(opt)

    def output_detail(self) -> None:
        """Output the detail to log file."""
        self.is_detail = True
        _, function = self._function.split()
        self._function = f'#P {function}\n'
        self._root_space = '  '

    def reset_variable(self) -> None:
        """Reset the variables of the gjf file."""
        self._root = []
        self._symbols = []
        self._coordinates = []
        self._link_dict = {}
        self._link_num = 0

    def set_charge_spin(self, charge: int, spin: int) -> None:
        """Set the charge and spin of the molecule.

        Parameters
        ----------
        charge : int
            Charge of the molecule.
        spin : int
            Spin of the molecule.
        """
        self._charge_spin = f'\n{charge} {spin}\n'

    def set_coordinates(self, coordinates: List[Tuple[float, float, float]]) -> None:
        """Set the coordinates of the molecule.

        Parameters
        ----------
        coordinates : List[Tuple[float, float, float]]
            Coordinates of the molecule.
        """
        self._coordinates.extend(coordinates)

    def set_function(self, function: str) -> None:
        """Set the function of the gjf file.

        Parameters
        ----------
        function : str
            Function of the gjf file.
        """
        detail = ''
        if self.is_detail:
            detail = 'P'

        self._function = f'#{detail} {function}\n'

    def set_resource(self, cpu_num: int, mem_num: int, mem_unit: str = 'GB') -> None:
        """Set the number of cpu and memory.

        Parameters
        ----------
        cpu_num : int
            Number of cpu.
        mem_num : int
            Number of memory.
        mem_unit : str
            Unit of memory, by default 'GB'.
        """
        self._cpu_num = f'%NProcShared={cpu_num}\n'
        self._mem_num = f'%Mem={mem_num}{mem_unit}\n'

    def set_symbols(self, symbols: List[str]) -> None:
        """Set the symbols of the molecule.

        Parameters
        ----------
        symbols : List[str]
            Symbols of the molecule.
        """
        self._symbols.extend(symbols)

    def set_title(self, title: str) -> None:
        """Set the title of the gjf file.

        Parameters
        ----------
        title : str
            Title of the gjf file.
        """
        self._title = f'\n{title}\n'

    def _check_charge_spin(self) -> None:
        """Check the charge and spin of the molecule."""
        elements = 0

        if not self._symbols:
            self._charge_spin = ''
            return

        for symbol in self._symbols:
            elements += self.ELEMENTS_NUM[symbol]

        if elements % 2 == 0:
            self._charge_spin = '\n0 1\n'
        else:
            self._charge_spin = '\n0 2\n'

    def _check_radical(self) -> bool:
        """Check if the molecule is radical.

        Returns
        -------
        bool
            True if the molecule is radical, False otherwise.

        Warnings
        --------
        This radical check is incomplete. If the number of unpaired electrons is even, it is not considered a radical.
        """
        # FIXME: This radical check is incomplete. If the number of unpaired electrons is even, it is not considered a radical.
        warnings.warn('This radical check is incomplete. If the number of unpaired electrons is even, it is not considered a radical.')
        elements = 0
        for symbol in self._symbols:
            elements += self.ELEMENTS_NUM[symbol]

        if elements % 2 == 0:
            return False
        else:
            return True
