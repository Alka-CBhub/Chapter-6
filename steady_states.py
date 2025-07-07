import numpy as np
from scipy.optimize import fsolve

def find_steady_states(system_func, num_vars, num_samples=500, tol=1e-8, domain=(0, 1)):
    """
    Finds steady-state solutions for an n-dimensional system using random initial guesses.

    Parameters:
    -----------
    system_func : callable
        A function f(X) that returns [dx1, dx2, ..., dxn].
    num_vars : int
        Number of variables in the system.
    num_samples : int, optional
        Number of random initial guesses. Default is 500.
    tol : float, optional
        Tolerance for considering two solutions the same. Default is 1e-8.
    domain : tuple of float, optional
        Bounds for each variable as (min, max). Default is (0, 1).

    Returns:
    --------
    roots : list of numpy arrays
        A list of unique steady-state vectors [x1_ss, ..., xn_ss] found in the domain.
    """
    roots = []
    rng = np.random.default_rng()
    guesses = rng.uniform(low=domain[0], high=domain[1], size=(num_samples, num_vars))

    for guess in guesses:
        sol, info, ier, _ = fsolve(system_func, guess, full_output=True)
        if ier != 1:
            continue

        # Check if solution lies within domain bounds
        if not np.all((domain[0] <= sol) & (sol <= domain[1])):
            continue

        # Check for duplicates
        if any(np.allclose(sol, existing, atol=tol) for existing in roots):
            continue

        roots.append(sol)

    return roots
