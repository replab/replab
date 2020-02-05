Error estimation in RepLAB
==========================

No iterative approximation error
--------------------------------

For computations that are supposed to be exact (as those obtained from the decomposition of finite groups), we assume (for now) that the error on the coefficients is within ``2*eps(coeff)``, where ``eps`` represents the distance to the neighboring floating point number.

A rough estimate of the Frobenius error on a resulting matrix ``X`` with that error model is ``norm(X, 'fro') * sqrt(size(X, 1) * size(X, 2)) * 5e-16``.

Presence of approximation error
-------------------------------

When using a fixed-point iteration to estimate projections, etc..., we fit an exponential function
plus a noise floor, and return the error estimated from the noise floor of the fit.
