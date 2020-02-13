Accompanying code for "Nearest \Omega-stable matrix via Riemannian optimization", V. Noferini and F. Poloni.

This code has been tested with Manopt 5.0.

Example:

    >> importmanopt
    Manopt was added to Matlab's path.
    >> rng(0)
    >> A = randn(5);
    >> B = nearest_real_eigenvalues(A);
    >> C = nearest_stable_real(A);     # nearest Hurwitz stable in R^(nxn)
    >> D = nearest_stable_complex(A);  # nearest Hurwitz stable in C^(nxn)

It may seem that the computed solutions do not respect the constraints, e.g.,

    >> eig(B)
    ans =
      -1.3732e+00 + 0.0000e+00i
       4.8213e-01 + 0.0000e+00i
       4.8200e-01 + 1.2719e-04i
       4.8200e-01 - 1.2719e-04i
       4.8187e-01 + 0.0000e+00i

does not have computed real eigenvalues. However, this is only due to round-off errors: the matrix B is at a distance comparable with machine precision from one with real eigenvalues. The fourth and fifth output value provide a Schur decomposition of the solution:

    >> [B, ~, ~, Q, T] = nearest_real_eigenvalues(A)
    B =
       9.7238e-01   8.7634e-01   1.4553e+00  -2.8240e-01  -2.6962e-01
       8.0089e-01  -1.0966e+00   2.6768e-01  -3.8187e-01  -1.1662e-01
      -3.9443e-01  -1.0929e+00  -7.2378e-01   1.1848e-01   3.1793e-01
       3.0884e-01  -7.9942e-01   1.3594e+00   2.8792e-01   1.1094e+00
      -8.2340e-01  -2.9726e+00  -1.6790e+00  -7.5791e-01   1.1149e+00
    Q =
      -3.8670e-01   6.5349e-01   1.9233e-01  -1.0284e-01  -6.1307e-01
       4.8406e-01   1.7892e-01   6.0460e-02  -8.5309e-01   4.7463e-02
       3.6523e-01  -4.9559e-01  -5.1955e-02   5.5984e-02  -7.8432e-01
      -3.8176e-01  -4.5324e-01   7.6290e-01  -2.5539e-01   3.9851e-02
       5.8053e-01   2.9984e-01   6.1208e-01   4.3966e-01   7.1710e-02
    T =
      -1.3732e+00   5.0898e-01  -3.9577e-01   2.5657e+00   2.0421e+00
                0   4.8200e-01  -7.0671e-01  -4.3808e-01  -5.9247e-01
                0            0   4.8200e-01   2.6800e+00  -2.7807e-01
                0            0            0   4.8200e-01   1.9519e+00
                0            0            0            0   4.8200e-01
    >> norm(B - Q*T*Q')
    ans =
         0

One can see in this example that T does indeed have real eigenvalues (including, in particular, one with multiplicity 4). It is consistent with eigenvalue perturbation theory that when one forms the product Q*T*Q' these eigenvalues are affected by perturbations of magnitude eps^(1/4) \approx 10^(-4).

The second and third output of the functions contain diagnostics in the same format as those produced by the code on http://bit.ly/NearestStableMatrix2 by N. Gillis (objective function and time after each iteration).
