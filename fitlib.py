import numpy as np
from numpy import linalg as LA
from inspect import signature
from scipy.optimize import curve_fit

class TFFit():
    def __init__(self, x, y, dy, fun, p0):
        self.x = x
        self.y = y
        self.dy = dy
        self.fun = fun
        self.npar = len(signature(fun).parameters)
        self.p0 = p0
        self.pars = {}

    def fit(self):
        popt, pcov = curve_fit(
                        self.fun,
                        self.x,
                        self.y,
                        sigma=self.dy,
                        p0=self.p0,
                        absolute_sigma=True,
                        maxfev=20000,
                        bounds=(0,np.inf)
                        )

        # Nominal best fit
        self.pars['best'] = popt
        
        # Diagonalize the covariance matrix
        pcov_w, pcov_v  = LA.eig(pcov)

        # Add the variations
        for i in range(len(pcov_w)):
            variation = np.sign(pcov_w[i])*np.sqrt(np.abs(pcov_w[i]))*pcov_v[:,i]
            self.pars[f'fit_{i}_dn'] = popt - variation
            self.pars[f'fit_{i}_up'] = popt + variation

    def envelope(self, x):
        vals = np.stack(self.evaluate_all(x).values())
        return np.min(vals, axis=0), np.max(vals, axis=0)


            
    def evaluate(self, x, variation='best'):
        return self.fun(x, *self.pars[variation])

    def evaluate_all(self, x):
        ret = {}
        for variation in self.pars.keys():
            ret[variation] = self.evaluate(x, variation)
        return ret
