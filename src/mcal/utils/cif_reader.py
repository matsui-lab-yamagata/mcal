"""CifReader beta (2025/10/30)"""
import os
import re
from itertools import product
from pathlib import Path
from typing import Dict, List, Literal, Tuple
import warnings

import numpy as np
import pandas as pd
from numpy.typing import NDArray


class CifReader:
    """CifReader class.

    This class is used to read cif file and extract crystal information.

    Raises
    ------
    ElementPropertiesIsNotDefinedError
        Raised when element properties is not defined.
    SymmetryIsNotDefinedError
        Raised when symmetry is not defined.
    ZValueIsNotMatchError
        Raised when z value is not match.
        The atomic bond detection may not be functioning correctly.
    """
    parent_dir = Path(os.path.abspath(__file__)).parent.parent
    ELEMENT_PROP = pd.read_csv(f'{parent_dir}/constants/element_properties.csv').dropna(axis=0)
    ATOMIC_WEIGHTS = ELEMENT_PROP[['symbol', 'weight']].set_index('symbol').to_dict()['weight']
    COVALENT_RADII = ELEMENT_PROP[['symbol', 'covalent_radius']].set_index('symbol').to_dict()['covalent_radius']

    def __init__(self, cif_path: str) -> None:
        """Initialize the CifReader class.

        Parameters
        ----------
        cif_path : str
            Path of cif file.
        """
        self.basename = None

        # Crystal properties
        self.cell_lengths = [None, None, None]
        self.cell_angles = [None, None, None]
        self.lattice = None
        self.symmetry_pos = []
        self.z_value = 0
        self._ref_z_value = 0

        # Molecule properties
        self.symbols = []
        self.symbols_label = []
        self.coordinates = []
        self.sym_symbols = []
        self.sym_coords = np.empty((0, 3))

        # Unique molecule
        self.unique_symbols = {}
        self.unique_coords = {}

        self._reader(cif_path)
        self._calc_lattice()
        self._operate_sym()
        self.sym_symbols, self.sym_coords = self.remove_duplicates(self.sym_symbols, self.sym_coords)
        self._make_adjacency_mat()
        self._split_mols()
        self._put_unit_cell()
        self.sym_symbols, self.sym_coords = self.remove_duplicates(self.sym_symbols, self.sym_coords)
        self._make_adjacency_mat()
        self._split_mols()
        self._calc_z_value()

        if self._ref_z_value != 0 and self.z_value != self._ref_z_value:
            raise ZValueIsNotMatchError('Z value is not match.')

    def _calc_lattice(self):
        """Calculate lattice."""
        a, b, c = self.cell_lengths
        alpha, beta, gamma = tuple(map(lambda x: np.radians(x), self.cell_angles))

        b_x = b * np.cos(gamma)
        b_y = b * np.sin(gamma)
        c_x = c * np.cos(beta)
        v = ((np.cos(alpha) - np.cos(beta) * np.cos(gamma))) / np.sin(gamma)
        c_y = c * v
        c_z = c * np.sqrt(1 - np.cos(beta)**2 - v**2)

        self.lattice = np.array((
            (a, 0, 0),
            (b_x, b_y, 0),
            (c_x, c_y, c_z),
        ))

    def _calc_z_value(self):
        """Calculate z value."""
        for atom_idx in self.bonded_atoms:
            cen_of_weight = self.calc_cen_of_weight(self.sym_coords[atom_idx])

            if self._is_in_unit_cell(cen_of_weight):
                self.unique_symbols[self.z_value] = self.sym_symbols[atom_idx]
                self.unique_coords[self.z_value] = self.sym_coords[atom_idx]
                self.z_value += 1

    def _is_in_unit_cell(self, cen_of_weight: NDArray[np.float64]) -> bool:
        """Determine if the center of weight is in the unit cell.

        Parameters
        ----------
        cen_of_weight : NDArray[np.float64]
            Center of weight.

        Returns
        -------
        bool
            True if the center of weight is in the unit cell.
        """
        if np.all(0 <= cen_of_weight) and np.all(cen_of_weight < 1):
            is_in_unit_cell = True
        else:
            is_in_unit_cell = False

        return is_in_unit_cell

    def _make_adjacency_mat(self):
        """Determine bonding and create the adjacency matrix."""
        num_atoms = len(self.sym_symbols)
        self.adjacency_mat = np.zeros((num_atoms, num_atoms), dtype=bool)

        self.cart_coords = np.dot(self.sym_coords, self.lattice)

        try:
            covalent_distance = np.array([self.COVALENT_RADII[symbol] for symbol in self.sym_symbols]) \
                + np.array([self.COVALENT_RADII[symbol] for symbol in self.sym_symbols])[:, np.newaxis]
        except KeyError:
            raise ElementPropertiesIsNotDefinedError('Element properties is not defined.')

        distance = np.linalg.norm(self.cart_coords[:, np.newaxis, :] - self.cart_coords[np.newaxis, :, :], axis=-1)
        self.adjacency_mat[(distance <= covalent_distance * 1.3) & (distance != 0)] = 1

    def _operate_sym(self) -> None:
        """Perform molecular symmetry operations."""

        def _extract_coord(coord: NDArray[np.float64], idx: int, is_minus: bool) -> NDArray[np.float64]:
            """Extract coordinates from the coordinate array.

            Parameters
            ----------
            coord : NDArray[np.float64]
                Coordinate array.
            idx : int
                Index of the coordinate to extract.
            is_minus : bool
                If True, the coordinate is extracted with a minus sign.

            Returns
            -------
            NDArray[np.float64]
                Extracted coordinate array.
            """
            if is_minus:
                return -coord[:, idx].copy()
            else:
                return coord[:, idx].copy()


        if len(self.symmetry_pos) == 0:
            raise SymmetryIsNotDefinedError('Symmetry is not defined.')

        self.sym_symbols = np.tile(self.symbols, len(self.symmetry_pos))

        idx_fil = ('x', 'y', 'z')

        for pos in self.symmetry_pos:
            moved_coord = np.zeros(self.coordinates.shape)

            sum_array = np.empty(0)

            for i, s, in enumerate(pos.split(',')):
                matches = re.findall(r'[0-9]/[0-9]', s)
                if matches:
                    fraction = eval(f'float({matches[0]})')
                else:
                    fraction = 0
                sum_array = np.append(sum_array, fraction)

                terms = [x.replace('+', '') for x in re.findall(r'\+?\-?[x-z]', s)]

                for term in terms:
                    is_minus = False

                    if '-' in term:
                        is_minus = True
                        term = term[-1]

                    moved_coord[:, i] += _extract_coord(
                        self.coordinates,
                        idx_fil.index(term),
                        is_minus
                    )
            moved_coord = sum_array + moved_coord
            self.sym_coords = np.append(self.sym_coords, moved_coord, axis=0)

    def _put_unit_cell(self) -> None:
        """Put molecules into unit cell."""
        for atom_idx in self.bonded_atoms:
            for i, c in enumerate(self.calc_cen_of_weight(self.sym_coords[atom_idx])):
                if 1 <= c:
                    change = -int(c)
                elif c < 0:
                    change = abs(int(c)) + 1
                else:
                    change = 0
                self.sym_coords[atom_idx, i] += change

    def _reader(self, cif_path: str) -> None:
        """Read cif file infomation.

        Parameters
        ----------
        cif_path : str
            Path of cif file.
        """
        # save index position
        counter = 0
        atom_data_index = {
            '_atom_site_label': None,
            '_atom_site_type_symbol': None,
            '_atom_site_fract_x': None,
            '_atom_site_fract_y': None,
            '_atom_site_fract_z': None,
        }
        symmetry_data_index = None

        is_read_atom = False
        is_read_sym = False

        with open(cif_path) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip()

                # remove blank characters
                if not line:
                    continue

                if line.startswith('data_'):
                    self.basename = '_'.join(line.split('_')[1:])

                # get unit cell information
                cell_params = ('_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma')
                if line.startswith(tuple(cell_params)):
                    value = float(re.sub(r'\(.*\)', '', line.split()[-1]))
                    if line.startswith('_cell_length'):
                        self.cell_lengths[cell_params.index(line.split()[0])] = value
                    else:
                        self.cell_angles[cell_params.index(line.split()[0])%3] = value
                elif line.startswith('_cell_formula_units_Z'):
                    self._ref_z_value = int(line.split()[-1])

                # get index position
                if 'loop_' == line:
                    counter = 0
                    is_read_atom = False
                    is_read_sym = False
                    continue
                elif '_' == line[0]:
                    if line in atom_data_index.keys():
                        atom_data_index[line] = counter
                        is_read_atom = True
                        is_read_sym = False
                    elif line in ('_symmetry_equiv_pos_as_xyz', '_space_group_symop_operation_xyz'):
                        symmetry_data_index = counter
                        is_read_atom = False
                        is_read_sym = True
                    else:
                        is_read_sym = False
                    counter += 1
                    continue
                elif ';' == line[0]:
                    is_read_atom = False
                    is_read_sym = False
                    continue

                if line[0] not in ('_', '#'):
                    # get symbol and fractional coordinates
                    if is_read_atom:
                        tmp_atom_data = line.split()
                        # remove disorder
                        if '?' not in tmp_atom_data[atom_data_index['_atom_site_label']]:
                            if atom_data_index['_atom_site_type_symbol'] is None:
                                symbol_label = tmp_atom_data[atom_data_index['_atom_site_label']]
                                symbol = symbol_label
                                for s in ['A', 'B', 'C']:
                                    symbol = symbol.replace(s, '')
                                symbol = re.sub(r'\d+', '', symbol)
                            else:
                                symbol_label = tmp_atom_data[atom_data_index['_atom_site_label']]
                                symbol = tmp_atom_data[atom_data_index['_atom_site_type_symbol']]
                            fract_x = tmp_atom_data[atom_data_index['_atom_site_fract_x']]
                            fract_y = tmp_atom_data[atom_data_index['_atom_site_fract_y']]
                            fract_z = tmp_atom_data[atom_data_index['_atom_site_fract_z']]
                            coord = [float(re.sub(r'\(.*\)', '', x)) for x in [fract_x, fract_y, fract_z]]
                            self.symbols.append(symbol)
                            self.symbols_label.append(symbol_label)
                            self.coordinates.append(coord)
                    # get symmetry operation information
                    elif is_read_sym:
                        if "'" in line:
                            line = list(map(lambda x: x.strip().replace(' ', ''), line.split("'")))
                            self.symmetry_pos.append(line[symmetry_data_index].lower())
                        else:
                            line = list(map(lambda x: x.strip().replace(' ', ''), line.split()))
                            self.symmetry_pos.append(line[symmetry_data_index].lower())

            self.symbols = np.array(self.symbols)
            self.coordinates = np.array(self.coordinates)

    def _search_connect_atoms(self, node: int, atoms: List[int], visited: NDArray[bool], num_atoms: int) -> None:
        """Find bonded atoms using depth-first search.

        Parameters
        ----------
        node : int
            Index of the atom.
        atoms : List[int]
            List of bonded atoms.
        visited : NDArray[bool]
            Array of visited atoms.
        num_atoms : int
            Number of atoms.
        """
        visited[node] = True
        atoms.append(node)
        for i in range(num_atoms):
            if self.adjacency_mat[node, i] and not visited[i]:
                self._search_connect_atoms(i, atoms, visited, num_atoms)

    def _split_mols(self) -> None:
        """Split molecules."""
        self.bonded_atoms = []
        num_atoms = len(self.sym_symbols)
        visited = np.zeros(num_atoms, dtype=bool)

        for i in range(num_atoms):
            if not visited[i]:
                atoms = []
                self._search_connect_atoms(i, atoms, visited, num_atoms)
                self.bonded_atoms.append(atoms)

        # get row corresponding to index
        self.sub_matrices = []
        for index_group in self.bonded_atoms:
            # get row corresponding to index
            index_group.sort()
            sub_matrix = self.adjacency_mat[np.ix_(index_group, index_group)]
            self.sub_matrices.append(sub_matrix)

    def calc_cen_of_weight(self, coordinates: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculate center of weight.

        Parameters
        ----------
        coordinates : NDArray[np.float64]
            Coordinates of monomolecular.

        Returns
        -------
        NDArray[np.float64]
            Center of weight.
        """
        cen_of_weight = np.average(coordinates, axis=0)

        return np.round(cen_of_weight, decimals=10)

    def convert_cart_to_frac(self, cart_coord: NDArray[np.float64]) -> NDArray[np.float64]:
        """Convert Cartesian coordinates to fractional coordinates.

        Parameters
        ----------
        cart_coord : NDArray[np.float64]
            Cartesian coordinates.

        Returns
        -------
        NDArray[np.float64]
            Fractional coordinates.
        """
        a, b, c = self.cell_lengths
        alpha, beta, gamma = tuple(map(lambda x: np.radians(x), self.cell_angles))
        b_x = -np.cos(gamma) / (a*np.sin(gamma))
        b_y = 1 / (b*np.sin(gamma))
        v = np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
        c_x = (np.cos(alpha)*np.cos(gamma) - np.cos(beta)) / (a*v*np.sin(gamma))
        c_y = (np.cos(beta)*np.cos(gamma) - np.cos(alpha)) / (b*v*np.sin(gamma))
        c_z = np.sin(gamma) / (c*v)

        vector = np.array((
                (1/a, 0, 0),
                (b_x, b_y, 0),
                (c_x, c_y, c_z),
            ))

        return np.dot(cart_coord, vector)

    def convert_frac_to_cart(self, frac_coord: NDArray[np.float64]) -> NDArray[np.float64]:
        """Convert fractional coordinates to Cartesian coordinates.

        Parameters
        ----------
        frac_coord : NDArray[np.float64]
            Fractional coordinates.

        Returns
        -------
        NDArray[np.float64]
            Cartesian coordinates.
        """
        return np.dot(frac_coord, self.lattice)

    def expand_mols(
        self,
        expand_range: int = 1
    ) -> Dict[Tuple[int, int, int], Dict[int, List[Tuple[str, NDArray[np.float64]]]]]:
        """Generate molecules around unique molecules.

        Parameters
        ----------
        expand_range : int
            The number of molecular cycles produced., by default 1

        Returns
        -------
        Dict[Tuple[int, int, int], Dict[int, List[Tuple[str, NDArray[np.float64]]]]]
            A nested dictionary containing the expanded molecular structure:

            - Outer key: Tuple[int, int, int]
                Represents the unit cell offset (i, j, k) relative to the origin unit cell.
                For example, (0, 0, 0) is the origin unit cell, (1, 0, 0) is one unit cell away in the a-direction, etc.

            - Inner key: int
                The index of the unique molecule within that unit cell.

            - Value: Tuple[List[str], NDArray[np.float64]]
                A list containing molecular information:
                - List[str]: Element symbols of the molecule
                - NDArray[np.float64]: Cartesian coordinates of the molecule (shape: (3, n))
        """
        expand_mols = {}
        combs = tuple(product(tuple(range(-expand_range, expand_range+1)), repeat=3))

        for comb in combs:
            for i, unique_coord in self.unique_coords.items():
                if i == 0:
                    expand_mols[comb] = {i: [self.unique_symbols[i], unique_coord + np.array(comb)]}
                else:
                    expand_mols[comb][i] = [self.unique_symbols[i], unique_coord + np.array(comb)]

        return expand_mols

    def export_unit_cell_file(self, file_path: str, format: Literal['mol', 'xyz'] = 'mol') -> None:
        """export unit cell file

        Parameters
        ----------
        file_path : str
            Path of the file to be saved.
        format : Literal['mol', 'xyz']
            Format of the file to be saved.
        """
        unit_cell_file = FileIO()
        for idx, symbols in self.unique_symbols.items():
            unit_cell_file.add_symbols(symbols)
            unit_cell_file.add_coordinates(self.convert_frac_to_cart(self.unique_coords[idx]))
            unit_cell_file.add_adjacency_mat(self.sub_matrices[idx])

        if format == 'mol' and unit_cell_file.atom_num > 999:
            format = 'xyz'
            file_path = file_path.replace('.mol', '.xyz')
            warnings.warn('The number of atoms is greater than 999. The file is saved as xyz format.')

        if format == 'mol':
            unit_cell_file.export_mol_file(file_path, header1=self.basename, header2="Generated by cif_reader.py")
        elif format == 'xyz':
            unit_cell_file.export_xyz_file(file_path, comment="Generated by cif_reader.py")

    def remove_duplicates(
        self,
        symbol: List[str],
        coordinate: NDArray[np.float64],
        tol: float = 1e-4,
    ) -> Tuple[List[str], NDArray[np.float64]]:
        """Remove duplicates from symbol and coordinate arrays based on coordinate with a given tolerance.

        Parameters
        ----------
        symbol : List[str]
            Symbols of molecules.
        coordinate : NDArray[np.float64]
            Coordinates of molecules.
        tol : float
            Tolerance for duplicate detection.

        Returns
        -------
        Tuple[List[str], NDArray[np.float64]]
            Symbols and coordinates of unique molecules.
        """
        distance_mat = ((coordinate[np.newaxis, :, :] - coordinate[:, np.newaxis, :]) ** 2).sum(axis=-1)
        dup = (distance_mat <= tol)
        dup = np.tril(dup, k=-1)
        unique_indices = ~dup.any(axis=-1)

        return symbol[unique_indices], coordinate[unique_indices]


class ElementPropertiesIsNotDefinedError(Exception):
    """Exception raised when element properties is not defined."""
    pass


class SymmetryIsNotDefinedError(Exception):
    """Exception raised when symmetry is not defined."""
    pass


class ZValueIsNotMatchError(Exception):
    """Exception raised when z value is not match."""
    pass


class FileIO:
    def __init__(self) -> None:
        self.atom_num = 0
        self.symbols_list = []
        self.coordinates_list = []
        self.adjacency_mat_list = []

    def add_adjacency_mat(self, adjacency_mat: NDArray[bool]) -> None:
        """add adjacency matrix

        Parameters
        ----------
        adjacency_mat : NDArray[bool]
            Adjacency matrix.
        """
        self.adjacency_mat_list.append(adjacency_mat)

    def add_coordinates(self, coordinates: NDArray[np.float64]) -> None:
        """add coordinates

        Parameters
        ----------
        coordinates : NDArray[np.float64]
            Coordinates.
        """
        self.coordinates_list.append(coordinates)

    def add_symbols(self, symbols: List[str]) -> None:
        """add symbols

        Parameters
        ----------
        symbols : List[str]
            Symbols.
        """
        self.atom_num += len(symbols)
        self.symbols_list.append(symbols)

    def export_mol_file(
        self,
        file_path: str,
        header1: str,
        header2: str,
    ) -> None:
        """export mol file

        Parameters
        ----------
        file_path : str
            Path of the file to be saved.
        header1 : str
            Header line 1.
        header2 : str
            Header line 2.
        """
        atom_lines = []
        bond_lines = []
        total_atoms = 0
        total_bonds = 0

        for i in range(len(self.symbols_list)):
            for s, (x, y, z) in zip(self.symbols_list[i], self.coordinates_list[i]):
                atom_lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {s:<2s}     0\n")

            for j in range(len(self.symbols_list[i])):
                for k in range(j):
                    if self.adjacency_mat_list[i][j, k] == 1:
                        bond_lines.append(f"{j+total_atoms+1:3d}{k+total_atoms+1:3d}  1\n")

            total_bonds += np.int32(np.sum(np.tril(self.adjacency_mat_list[i], k=-1)))
            total_atoms += len(self.symbols_list[i])

        with open(file_path, 'w') as f:
            f.write(f"{header1}\n")
            f.write(f"{header2}\n")
            f.write("\n")

            f.write(f"{total_atoms:3d}{total_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n")

            f.writelines(atom_lines)
            f.writelines(bond_lines)

            f.write("M  END\n")
            f.write("$$$$\n")

    def export_xyz_file(
        self,
        file_path: str,
        comment: str,
    ) -> None:
        """export xyz file

        Parameters
        ----------
        file_path : str
            Path of the file to be saved.
        comment : str
            Comment.
        """
        xyz_file_lines = []
        total_atoms = 0
        for i in range(len(self.symbols_list)):
            for s, (x, y, z) in zip(self.symbols_list[i], self.coordinates_list[i]):
                xyz_file_lines.append(f"{s:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

            total_atoms += len(self.symbols_list[i])

        with open(file_path, 'w') as f:
            f.write(f"{total_atoms:3d}\n")
            f.write(f"{comment}\n")
            f.writelines(xyz_file_lines)
