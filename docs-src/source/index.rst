==================================
mcal Documentation
==================================

.. image:: https://img.shields.io/badge/python-3.9%20or%20newer-blue
   :target: https://www.python.org
   :alt: Python

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

Overview
========

``mcal`` is a tool for calculating mobility tensors of organic semiconductors.
It calculates transfer integrals and reorganization energy from crystal structures (CIF files),
and determines mobility tensors considering anisotropy and path continuity.

Features
--------

* Read crystal structures from CIF files
* Quantum chemical calculations using Gaussian 09/16 or PySCF
* Calculation of transfer integrals
* Calculation of reorganization energy
* Computation of mobility tensors and eigenvalues

Requirements
============

* Python 3.9 or newer
* NumPy
* Pandas
* Matplotlib
* yu-tcal==4.0.3

Quantum Chemistry Calculation Tools
------------------------------------

At least one of the following is required:

* Gaussian 09 or 16
* PySCF (macOS / Linux / WSL2(Windows Subsystem for Linux))
* GPU4PySCF (macOS / Linux / WSL2(Windows Subsystem for Linux))

.. important::
   When using Gaussian, the path to Gaussian must be set in your environment.
   PySCF is supported on macOS / Linux. Windows users must use WSL2.

Installation
============

Using Gaussian 09 or 16 (without PySCF)
----------------------------------------

.. code-block:: bash

   pip install yu-mcal

Using PySCF (CPU only, macOS / Linux / WSL2)
---------------------------------------------

.. code-block:: bash

   pip install "yu-mcal[pyscf]"

Using GPU acceleration with PySCF (macOS / Linux / WSL2)
---------------------------------------------------------

1. Check your installed CUDA Toolkit version:

.. code-block:: bash

   nvcc --version

2. Install mcal with GPU acceleration:

If your CUDA Toolkit version is 12.x:

.. code-block:: bash

   pip install "yu-mcal[gpu4pyscf-cuda12]"

If your CUDA Toolkit version is 11.x:

.. code-block:: bash

   pip install "yu-mcal[gpu4pyscf-cuda11]"


Verify Installation
-------------------

After installation, you can verify by running:

.. code-block:: bash

   mcal --help

mcal Usage Manual
=================

Basic Usage
-----------

.. code-block:: bash

   mcal <cif_filename or pkl_filename> <osc_type> [options]

Required Arguments
~~~~~~~~~~~~~~~~~~

* ``cif_filename``: Path to the CIF file
* ``pkl_filename``: Path to the pickle file
* ``osc_type``: Organic semiconductor type

  * ``p``: p-type semiconductor (uses HOMO level)
  * ``n``: n-type semiconductor (uses LUMO level)

Basic Examples
~~~~~~~~~~~~~~

.. code-block:: bash

   # Calculate as p-type semiconductor
   mcal xxx.cif p

   # Calculate as n-type semiconductor
   mcal xxx.cif n

Options
-------

Calculation Settings
~~~~~~~~~~~~~~~~~~~~

``-M, --method <method>``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Specify the calculation method (used in both Gaussian and PySCF calculations).

* **Default**: ``B3LYP/6-31G(d,p)``
* **Example**: ``mcal xxx.cif p -M "B3LYP/6-31G(d)"``

``-c, --cpu <number>``
^^^^^^^^^^^^^^^^^^^^^^

Specify the number of CPUs to use.

* **Default**: ``4``
* **Example**: ``mcal xxx.cif p -c 8``

``-m, --mem <memory>``
^^^^^^^^^^^^^^^^^^^^^^

Specify the amount of memory in GB.

* **Default**: ``10``
* **Example**: ``mcal xxx.cif p -m 16``

``-g, --g09``
^^^^^^^^^^^^^

Use Gaussian 09 (default is Gaussian 16).

* **Example**: ``mcal xxx.cif p -g``

PySCF Settings
~~~~~~~~~~~~~~

``--pyscf``
^^^^^^^^^^^^

Use PySCF instead of Gaussian for all calculations. Requires ``yu-mcal[pyscf]``.

* **Example**: ``mcal xxx.cif p --pyscf``

``--gpu4pyscf``
^^^^^^^^^^^^^^^^

Use GPU acceleration via gpu4pyscf. Automatically enables PySCF mode (no need to specify ``--pyscf``).
Requires ``yu-mcal[gpu4pyscf-cuda11]`` or ``yu-mcal[gpu4pyscf-cuda12]``.

* **Example**: ``mcal xxx.cif p --gpu4pyscf``

``--cart``
^^^^^^^^^^^

Use Cartesian basis functions instead of spherical harmonics (PySCF only).

* **Example**: ``mcal xxx.cif p --pyscf --cart``

Calculation Control
~~~~~~~~~~~~~~~~~~~

``-r, --read``
^^^^^^^^^^^^^^

Read results from existing files without executing calculations.
With Gaussian, reads from log files; with PySCF, reads from checkpoint (``.chk``) files.

* **Example**: ``mcal xxx.cif p -r``

``-rp, --read_pickle``
^^^^^^^^^^^^^^^^^^^^^^

Read results from existing pickle file without executing calculations.

* **Example**: ``mcal xxx_result.pkl p -rp``

``--resume``
^^^^^^^^^^^^

Resume calculation using existing results.
With Gaussian, checks log file termination; with PySCF, checks for existing checkpoint (``.chk``) files.

* **Example**: ``mcal xxx.cif p --resume``

``--fullcal``
^^^^^^^^^^^^^

Disable all speedup processing and calculate transfer integrals for all pairs from scratch.
The following two optimizations are disabled:

1. **Pair screening**: pairs are normally skipped based on moment of inertia and center-of-mass
   distance; ``--fullcal`` disables this screening.
2. **Monomer caching**: monomer SCF calculations for the same molecule type are normally skipped
   by reusing a previously computed result file; ``--fullcal`` forces all monomer calculations
   to be performed from scratch.

* **Example**: ``mcal xxx.cif p --fullcal``

``--no-monomer-cache``
^^^^^^^^^^^^^^^^^^^^^

Disable only monomer caching. Pair screening remains active.
All monomer SCF calculations are performed from scratch instead of reusing
previously computed result files.
When performing detailed transfer integral analysis using
`tcal <https://github.com/matsui-lab-yamagata/tcal>`_, it is recommended to use this option.

* **Example**: ``mcal xxx.cif p --no-monomer-cache``

``--cellsize <number>``
^^^^^^^^^^^^^^^^^^^^^^^

Specify the number of unit cells to expand in each direction around the central unit cell
for transfer integral calculations.

* **Default**: ``2`` (creates 5x5x5 supercell)
* **Examples**:

  * ``mcal xxx.cif p --cellsize 1`` (creates 3x3x3 supercell)
  * ``mcal xxx.cif p --cellsize 3`` (creates 7x7x7 supercell)

Output Settings
~~~~~~~~~~~~~~~

``-p, --pickle``
^^^^^^^^^^^^^^^^

Save calculation results to a pickle file.

* **Example**: ``mcal xxx.cif p -p``

``--plot-plane <plane>``
^^^^^^^^^^^^^^^^^^^^^^^^^

Plot mobility tensor as a 2D polar plot on specified crystallographic plane.

* **Available planes**: ``ab``, ``ac``, ``ba``, ``bc``, ``ca``, ``cb``
* **Default**: None (no plot generated)
* **Examples**:

  * ``mcal xxx.cif p --plot-plane ab`` (plot on ab-plane)
  * ``mcal xxx.cif p --plot-plane bc`` (plot on bc-plane)

Practical Usage Examples
-------------------------

Basic Calculations
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Calculate mobility of p-type xxx
   mcal xxx.cif p

   # Use 8 CPUs and 16GB memory
   mcal xxx.cif p -c 8 -m 16

High-Precision Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Calculate transfer integrals for all pairs (high precision, time-consuming)
   mcal xxx.cif p --fullcal

   # Use larger supercell to widen transfer integral calculation range
   mcal xxx.cif p --cellsize 3

   # Use different basis set
   mcal xxx.cif p -M "B3LYP/6-311G(d,p)"

PySCF Calculations
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Calculate using PySCF (CPU)
   mcal xxx.cif p --pyscf

   # Calculate using PySCF with GPU acceleration (no --pyscf needed)
   mcal xxx.cif p --gpu4pyscf

   # Use 8 CPUs and 16GB memory with PySCF
   mcal xxx.cif p --pyscf -c 8 -m 16

   # Resume interrupted PySCF calculation
   mcal xxx.cif p --pyscf --resume

   # Read from existing PySCF checkpoint files
   mcal xxx.cif p --pyscf -r

Reusing Results
~~~~~~~~~~~~~~~

.. code-block:: bash

   # Read from existing calculation results
   mcal xxx.cif p -r

   # Read from existing pickle file
   mcal xxx_result.pkl p -rp

   # Resume interrupted calculation
   mcal xxx.cif p --resume

   # Save results to pickle file
   mcal xxx.cif p -p

Output
------

Standard Output
~~~~~~~~~~~~~~~

* Reorganization energy
* Transfer integrals for each pair
* Diffusion coefficient tensor
* Mobility tensor
* Eigenvalues and eigenvectors of mobility

Generated Files
~~~~~~~~~~~~~~~

Reorganization Energy Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following files are generated during reorganization energy calculation
(``c`` = cation for p-type, ``a`` = anion for n-type):

**Gaussian:**

* ``xxx_opt_n.gjf`` / ``xxx_opt_n.log`` — geometry optimization of neutral molecule
* ``xxx_c.gjf`` / ``xxx_c.log`` (or ``xxx_a``) — SP energy of ion at neutral geometry
* ``xxx_opt_c.gjf`` / ``xxx_opt_c.log`` (or ``xxx_opt_a``) — geometry optimization of ion
* ``xxx_n.gjf`` / ``xxx_n.log`` — SP energy of neutral at ion geometry

**PySCF:**

* ``xxx_opt_n.xyz`` / ``xxx_opt_n.chk`` — geometry optimization of neutral molecule
* ``xxx_c.chk`` (or ``xxx_a.chk``) — SP energy of ion at neutral geometry
* ``xxx_opt_c.xyz`` / ``xxx_opt_c.chk`` (or ``xxx_opt_a``) — geometry optimization of ion
* ``xxx_n.chk`` — SP energy of neutral at ion geometry

Transfer Integral Files
^^^^^^^^^^^^^^^^^^^^^^^

mcal generates calculation files named using the ``(s_t_i_j_k)`` notation:

.. list-table::
   :header-rows: 1

   * - Symbol
     - Meaning
   * - ``s``
     - Molecule index in the reference unit cell (0,0,0)
   * - ``t``
     - Molecule index in the neighboring unit cell
   * - ``i``
     - Translation index along the **a**-axis
   * - ``j``
     - Translation index along the **b**-axis
   * - ``k``
     - Translation index along the **c**-axis

**Example:** ``xxx-(0_0_1_0_0)`` represents the transfer integral between the 0th molecule
in the (0,0,0) cell and the 0th molecule in the (1,0,0) cell.

**Gaussian:**

* ``xxx-(s_t_i_j_k).gjf`` / ``xxx-(s_t_i_j_k).log`` — dimer
* ``xxx-(s_t_i_j_k)_m1.gjf`` / ``xxx-(s_t_i_j_k)_m1.log`` — monomer 1
* ``xxx-(s_t_i_j_k)_m2.gjf`` / ``xxx-(s_t_i_j_k)_m2.log`` — monomer 2

**PySCF:**

* ``xxx-(s_t_i_j_k).xyz`` / ``xxx-(s_t_i_j_k).chk`` — dimer
* ``xxx-(s_t_i_j_k)_m1.chk`` — monomer 1
* ``xxx-(s_t_i_j_k)_m2.chk`` — monomer 2

Notes
-----

.. note::

   1. **Calculation Time**: Calculation time varies significantly depending on the number of
      molecules and cell size. By default, two speedup mechanisms are enabled: pair pre-screening
      (skipping pairs unlikely to have significant transfer integrals) and monomer caching
      (reusing the isolated-molecule SCF result for molecule types already computed).
      Use ``--fullcal`` to disable both.
   2. **Memory Usage**: Ensure sufficient memory for large systems
   3. **Gaussian Installation**: Gaussian 09 or Gaussian 16 is required
   4. **Dependencies**: Make sure all required Python libraries are installed

Troubleshooting
---------------

If calculation stops midway
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Resume with --resume option
   mcal xxx.cif p --resume

Memory shortage error
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Increase memory amount
   mcal xxx.cif p -m 32

To reduce calculation time
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Enable speedup processing (default)
   mcal xxx.cif p

   # Use smaller supercell for faster calculation
   mcal xxx.cif p --cellsize 1

   # Increase number of CPUs
   mcal xxx.cif p -c 16

API Reference
=============

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   modules

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Authors
=======

`Matsui Laboratory, Research Center for Organic Electronics (ROEL), Yamagata University <https://matsui-lab.yz.yamagata-u.ac.jp/index-e.html>`_

* Hiroyuki Matsui
* Koki Ozawa

Email: h-matsui[at]yz.yamagata-u.ac.jp (Please replace ``[at]`` with ``@``)

License
=======

This project is released under the MIT License.

For more details, see the `LICENSE file on GitHub <https://github.com/matsui-lab-yamagata/mcal/blob/main/LICENSE>`_.

Acknowledgements
================

This work was supported by JSPS Grant-in-Aid for JSPS Fellows Grant Number JP25KJ0647.
