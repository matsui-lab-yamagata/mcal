"""hopping_mobility_model.py (2025/10/06)"""
import math
import random
from typing import List, Tuple

import numpy as np
from numpy.typing import NDArray


const_kb = 1.380649e-23  # Boltzmann constant [J/K]
const_e = 1.60217663e-19  # Elementary charge [C(=J/eV)]
const_hbar = 6.62607015e-34 / (2 * math.pi)  # Dirac constant [Js]


def demo():
    # 一次元系で、粒子は1 sごとに0.01の確率で右へ1 m、0.01の確率で左へ1 m移動し、0.98の確率でその場に留まる。
    print("\nOne-dimensional system")
    lattice = np.array(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
    hop = ((0, 0, 1, 0, 0, 0.01),)
    D = diffusion_coefficient_tensor(lattice, hop)
    print("Diffusion coefficient tensor (analytical):")
    for d in D:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")
    print("Diffusion coefficient tensor (ODE):")
    D_ode = diffusion_coefficient_tensor_ODE(lattice, hop)
    for d in D_ode:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")
    print("Diffusion coefficient tensor (MC):")
    D_mc = diffusion_coefficient_tensor_MC(lattice, hop)
    for d in D_mc:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")

    # 一次元系で、粒子は偶数サイトにいる時は1sごとに0.02の確率で右へ1m、0.01の確率で左へ1m移動し、0.97の確率でその場に留まる。
    # 奇数サイトにいる時は1sごとに0.01の確率で右へ1m、0.02の確率で左へ1m移動し、0.97の確率でその場に留まる。
    print("\nOne-dimensional dimer system")
    lattice = np.array(((2.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
    hop = (
        (0, 1, 0, 0, 0, 0.02),
        (0, 1, -1, 0, 0, 0.01),
    )
    D = diffusion_coefficient_tensor(lattice, hop)
    print("Diffusion coefficient tensor (analytical):")
    for d in D:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")
    print("Diffusion coefficient tensor (ODE):")
    D_ode = diffusion_coefficient_tensor_ODE(lattice, hop)
    for d in D_ode:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")
    print("Diffusion coefficient tensor (MC):")
    D_mc = diffusion_coefficient_tensor_MC(lattice, hop)
    for d in D_mc:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")

    # ランダムな条件で検証
    print("\nRandom system")
    lattice = np.random.random((3, 3)) * 2.0
    lattice[0, 1:] = 0.0
    lattice[1, 2:] = 0.0
    hop = []
    for _ in range(6):
        s = random.randint(0, 1)
        t = random.randint(0, 1)
        i = random.randint(-1, 1)
        j = random.randint(-1, 1)
        k = random.randint(-1, 1)
        p = random.random() * 0.02
        hop.append((s, t, i, j, k, p))
    D = diffusion_coefficient_tensor(lattice, hop)
    print("Diffusion coefficient tensor (analytical):")
    for d in D:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")
    print("Diffusion coefficient tensor (ODE):")
    D_ode = diffusion_coefficient_tensor_ODE(lattice, hop)
    for d in D_ode:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")
    print("Diffusion coefficient tensor (MC):")
    D_mc = diffusion_coefficient_tensor_MC(lattice, hop)
    for d in D_mc:
        print(f"{d[0]:9.6f} {d[1]:9.6f} {d[2]:9.6f}")


def cal_pinv(array: NDArray[np.float64], rcond: float = 1e-9) -> NDArray[np.float64]:
    """Calculate pseudo-inverse matrix using eigenvalue decomposition

    Parameters
    ----------
    array : NDArray[np.float64]
        Input matrix
    rcond : float, optional
        Cutoff for small singular values, by default 1e-9

    Returns
    -------
    NDArray[np.float64]
        Pseudo-inverse matrix

    Raises
    ------
    ValueError
        The last eigenvalue is not zero.
    ValueError
        All eigenvalues except the last one should be negative.
    """
    eigvals, eigvecs = np.linalg.eigh(array)

    # Calculate pseudo-inverse matrix using eigenvalue decomposition
    inveigvals = np.zeros_like(eigvals)
    if abs(eigvals[-1] / eigvals[0]) > rcond:
        raise ValueError(f"The last eigenvalue is not zero, which is unexpected for this test case. {eigvals}")
    if any(eigvals[0:-1] > 0):
        raise ValueError(f"All eigenvalues except the last one should be negative, which is unexpected for this test case. {eigvals}")

    inveigvals[0:-1] = 1.0 / eigvals[0:-1]
    inveigvals[-1] = 0.0
    pinv = eigvecs @ np.diag(inveigvals) @ eigvecs.T

    return pinv


def marcus_rate(transfer: float, reorganization: float, T: float = 300.0) -> float:
    """Calculate hopping rate (1/s) from transfer integral (eV) and reorganization energy (eV)

    Parameters
    ----------
    transfer : float
        Transfer integral [eV]
    reorganization : float
        Reorganization energy [eV]
    T : float
        Temperature [K], by default 300.0

    Returns
    -------
    float
        Hopping rate [1/s]
    """
    kbT = const_kb * T
    return (
        (transfer * const_e) ** 2
        / const_hbar
        * math.sqrt(math.pi / (reorganization * const_e * kbT))
        * math.exp(-reorganization * const_e / (4 * kbT))
    )


def mobility_tensor(D: NDArray[np.float64], T: float = 300.0) -> NDArray[np.float64]:
    """Calculate mobility tensor from diffusion coefficient tensor

    Parameters
    ----------
    D : 3x3 numpy.array
        Diffusion coefficient tensor
    T : float
        Temperature [K], by default 300.0

    Returns
    -------
    3x3 numpy.array
        Mobility tensor
    """
    return D * const_e / (const_kb * T)


def diffusion_coefficient_tensor(
    lattice: NDArray[np.float64],
    hop: List[Tuple[int, int, int, int, int, float]]
) -> NDArray[np.float64]:
    """Calculate diffusion coefficient tensor from hopping rate

    Parameters
    ----------
    lattice : 3x3 numpy.array
        lattice[0,:] is a-axis vector, lattice[1,:] b-axis vector, lattice[2,:] c-axis vector
    hop : list of (int, int, int, int, int, float) tuple.
        (s, t, i, j, k, p) means that the hopping rate from s-th molecule in (0, 0, 0) cell to t-th molecule in (i, j, k) cell is p.

    Returns
    -------
    3x3 numpy.array
        Diffusion coefficient tensor
    """

    # Standardize hop list
    hop = _standardize_hop_list(hop)

    # Number of molecules in the unit cell
    n = len(set([h[0] for h in hop]) | set([h[1] for h in hop]))

    # Prepare arrays
    D = np.zeros((3, 3))
    B = np.zeros((n, n))
    C = np.zeros((n, 3))

    for s, t, i, j, k, p in hop:
        vec = np.array((i, j, k)) @ lattice
        D[:, :] += p * np.outer(vec, vec) * 2  # Consider hopping in both directions
        B[s, t] += p
        B[t, s] += p  # Consider hopping in both directions
        B[s, s] -= p
        B[t, t] -= p
        C[s, :] += p * vec
        C[t, :] -= p * vec  # Consider hopping in both directions

    # For n = 1 case, skip C.T @ B_pinv @ C term as it equals zero
    if n > 1:
        B_pinv = cal_pinv(B)
        D = (D / 2 + C.T @ B_pinv @ C) / n
    else:
        D = D / 2

    # Check computational errors
    threshold = np.max(abs(D)) * 1e-6
    D_diff = abs(D - D.T)
    if np.any(D_diff > threshold):
        raise ValueError(f"Diffusion coefficient tensor D should be symmetric: {D}")

    # Make symmetric matrix considering computational errors
    D = (D + D.T) / 2

    return D


def diffusion_coefficient_tensor_ODE(
    lattice: NDArray[np.float64],
    hop: List[Tuple[int, int, int, int, int, float]],
    max_steps: int = 200,
    size: int = 40,
    max_rate: float = 0.05
) -> NDArray[np.float64]:
    """Calculate diffusion coefficient tensor from numerical solution of Ordinary Differential Equation (ODE)

    Parameters
    ----------
    lattice : 3x3 numpy.array
        lattice[0,:] is a-axis vector, lattice[1,:] b-axis vector, lattice[2,:] c-axis vector
    hop : list of (int, int, int, int, int, float) tuple.
        (s, t, i, j, k, p) means that the hopping rate from s-th molecule in (0, 0, 0) cell to t-th molecule in (i, j, k) cell is p.
    max_steps : int
        Maximum number of steps
    size : int
        Size of the simulation box
    max_rate : float
        Maximum rate of hopping

    Returns
    -------
    3x3 numpy.array
        Diffusion coefficient tensor
    """
    # Standardize hop list
    hop = _standardize_hop_list(hop)

    # Number of molecules in the unit cell
    n = len(set([h[0] for h in hop]) | set([h[1] for h in hop]))

    prob = np.zeros((2 * size, 2 * size, 2 * size, n))
    prob[size, size, size, :] = 1.0 / n
    pre_prob = np.zeros_like(prob)

    dt = max_rate / max(h[5] for h in hop)  # Time step
    for _ in range(max_steps):
        pre_prob[:, :, :, :] = prob
        for s, t, i, j, k, p in hop:
            prob[:, :, :, t] += np.roll(pre_prob[:, :, :, s], (i, j, k), axis=(0, 1, 2)) * p * dt
            prob[:, :, :, s] += np.roll(pre_prob[:, :, :, t], (-i, -j, -k), axis=(0, 1, 2)) * p * dt
            prob[:, :, :, s] -= pre_prob[:, :, :, s] * p * dt
            prob[:, :, :, t] -= pre_prob[:, :, :, t] * p * dt

    # Check the sum of probabilities
    total_prob = np.sum(prob)
    assert np.isclose(total_prob, 1.0), f"Total probability is not 1: {total_prob}"

    # Average of outer products of positions
    avg_outer_product = np.zeros((3, 3))
    for i, j, k, l in np.ndindex(2 * size, 2 * size, 2 * size, n):
        vec = np.array((i - size, j - size, k - size)) @ lattice
        avg_outer_product += prob[i, j, k, l] * np.outer(vec, vec)

    D = avg_outer_product / (2 * max_steps * dt)
    return D


def diffusion_coefficient_tensor_MC(
    lattice: NDArray[np.float64],
    hop: List[Tuple[int, int, int, int, int, float]],
    steps: int = 100,
    particles: int = 10000
) -> NDArray[np.float64]:
    """Calculate diffusion coefficient tensor from Monte Carlo simulation using Gillespie algorithm.

    Parameters
    ----------
    lattice : 3x3 numpy.array
        lattice[0,:] is a-axis vector, lattice[1,:] b-axis vector, lattice[2,:] c-axis vector
    hop : list of (int, int, int, int, int, float) tuple.
        (s, t, i, j, k, p) means that the hopping rate from s-th molecule in (0, 0, 0) cell to t-th molecule in (i, j, k) cell is p.
    steps : int
        Number of steps
    particles : int
        Number of particles

    Returns
    -------
    3x3 numpy.array
        Diffusion coefficient tensor
    """
    # Standardize hop list
    hop = _standardize_hop_list(hop)

    # Number of molecules in the unit cell
    n = len(set([h[0] for h in hop]) | set([h[1] for h in hop]))

    paths = [[] for _ in range(n)]
    probs = [[] for _ in range(n)]
    total_rates = [0] * n
    for s, t, i, j, k, p in hop:
        paths[s].append((t, i, j, k))
        paths[t].append((s, -i, -j, -k))
        probs[s].append(p)
        probs[t].append(p)
        total_rates[s] += p
        total_rates[t] += p
    max_time = steps / np.mean(total_rates)

    # Simulation
    sum_outer_product = np.zeros((3, 3))
    for _ in range(particles):
        xyz = np.zeros(3, dtype=int)
        mol = random.choice(range(n))  # ランダムに初期位置を選ぶ
        t = 0.0
        while t < max_time:
            t += -math.log(1.0 - random.random()) / total_rates[mol]
            path = random.choices(paths[mol], weights=probs[mol])[0]
            xyz += path[1], path[2], path[3]
            mol = path[0]
        xyz = xyz @ lattice
        sum_outer_product += np.outer(xyz, xyz)

    # Calculate diffusion coefficient
    D = sum_outer_product / (particles * 2 * max_time)
    return D


def print_tensor(tensor: NDArray[np.float64], msg: str = 'Mobility tensor'):
    print('-' * (len(msg)+2))
    print(f' {msg} ')
    print('-' * (len(msg)+2))
    if tensor.shape == (3, ):
        print(f"{tensor[0]:12.6g} {tensor[1]:12.6g} {tensor[2]:12.6g}")
    elif tensor.shape == (3, 3):
        for a in tensor:
            print(f"{a[0]:12.6g} {a[1]:12.6g} {a[2]:12.6g}")
    print()


def _standardize_hop_list(hop: List[Tuple[int, int, int, int, int, float]]) -> List[Tuple[int, int, int, int, int, float]]:
    """
    Standardize the hop list by ensuring that s <= t.
    If s == t, ensure that the first non-zero component of (i, j, k) is positive.

    Parameters
    ----------
    hop: list of (int, int, int, int, int, float) tuples.
        (s, t, i, j, k, p) means that the hopping rate from s-th molecule in (0, 0, 0) cell to t-th molecule in (i, j, k) cell is p.

    Returns
    -------
    list of (int, int, int, int, int, float) tuple.
        List of standardized hopping rate tuples.
    """
    hop = list(hop)
    standardized_hop = []
    standardized_hop_keys = set()
    for s, t, i, j, k, p in hop:
        if s > t:
            s, t, i, j, k, p = t, s, -i, -j, -k, p
        elif s == t:
            if (i, j, k) == (0, 0, 0):
                continue
            elif i < 0 or (i == 0 and (j < 0 or (j == 0 and k < 0))):
                s, t, i, j, k, p = t, s, -i, -j, -k, p

        if (s, t, i, j, k) not in standardized_hop_keys:
            standardized_hop_keys.add((s, t, i, j, k))
            standardized_hop.append((s, t, i, j, k, p))

    return standardized_hop


if __name__ == "__main__":
    demo()
