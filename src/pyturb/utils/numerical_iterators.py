"""

MRodriguez 2020

"""
import warnings
import numpy as np

root_tolerance = 1e-8
root_niter_max = 1000


def variable_step_roots(x0, func, dxmax=1, verbosity=False, root_niter_max_=root_niter_max):

    itera = True
    dx = dxmax
    factor = 1

    residual = func(x0)
    xi = x0 + dx
    ii = 1

    while itera:
        residual0 = residual
        residual = func(xi)

        if verbosity:
            print(ii, xi, dx, residual, residual0)

        if np.abs(dx)<root_tolerance:
            itera=False

        elif np.sign(residual) == np.sign(residual0):
            if np.abs(residual)<np.abs(residual0):
                x0 = xi
                xi += dx

            else:
                dx=-dx
                x0 = xi
                xi += dx

        else:
            factor *= 0.5
            dx = -dx*factor
            xi += dx

        if ii==root_niter_max_:
            itera = False
            warnings.warn("Maximum iterations allowed exceeded: {0}".format(root_niter_max_),UserWarning)

        else:
            ii += 1

            
    return xi, dx, residual
