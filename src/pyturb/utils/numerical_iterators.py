"""

MRodriguez 2020

"""
import warnings
import numpy as np



def variable_step_roots(x0, func, dxmax=1, verbosity=False, root_niter_max_=1000, root_tolerance=1e-6):
    """
    """
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


def two_equations_system_guess(x0A, eqA, eqB, verbosity=False, niter_max=100, tolerance=1e-8):
    """
    """


    iterate = True
    it_counter = 0

    residual0 = 0
    
    while iterate:
        it_counter += 1

        xiB_guess = eqA(x0A)
        xiA = eqB(xiB_guess)

        residual = ((xiA - x0A)/x0A)

        if verbosity:
            print('---',it_counter, x0A, xiA, residual)

        if it_counter == 1:
            residual0 = residual
            x0A = xiA
            pass
        elif np.abs(residual)<tolerance:
            iterate = False
            
        elif np.sign(residual)==np.sign(residual0):
            if np.abs(residual)>np.abs(residual0):
                xiA = x0A + (xiA-x0A)/2
                print('   -B-', residual, residual0, xiA)
            else:
                x0A = xiA
                print('   -A-', residual, residual0, xiA)
                residual0 = residual
            
        else:
            xiA = x0A + (xiA-x0A)/2
            print('   -C-', residual, residual0, xiA)

        if it_counter==niter_max:
            iterate = False
        


    return xiA
