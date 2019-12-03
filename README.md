# AnySOS
This is a MATLAB implementation of the code from the paper "AnySOS: An anytime algorithm for semidefinite programming" by Driggs and Fawzi (CDC, 2019). AnySOS is suited for solving large-scale semidefinite programs where suboptimal solutions are acceptable, but feasibility in the solution is required. Example problems of this type include sum-of-squares problems arising in control applications.

If you wish to cite AnySOS please use the following:
```
@inproceedings{anysos2019,
    author       = {D. Driggs and H. Fawzi},
    title        = {{AnySOS}: {A}n anytime algorithm for semidefinite programming},
    journal      = {Conference on Decision and Control},
    month        = {December},
    year         = {2019}
}
@misc{anysos,
    author       = {D. Driggs and H. Fawzi},
    title        = {{AnySOS}, version 1.0.0},
    howpublished = {\url{https://github.com/derekdriggs/AnySOS}},
    month        = dec,
    year         = 2019
}
```
AnySOS solves problems of the form
```
minimize        c'x
subject to      Ax = b
                x in K
```
where the cone `K` can be any Cartesian product of the following:
+ free cone `{x | x in R}`)
+ positive orthant `{x | x >= 0}`
+ positive semidefinite cone `{ X | min(eig(X)) >= 0, X = X^T }`

The columns of `A` must correspond to variables in the above cones in the given order (i.e. the first columns of `A` correspond to the free variables, then the non-negative variables, and finally the semidefinite variables.) For semidefinite variables, AnySOS uses the convention that the matrix variables `X` are vectorized by scaling the off-diagonal entries by `sqrt(2)` and stacking the lower triangle using column-major indexing. The function `my_svec(X)` performs this operation, and `my_smat(x)` is its inverse.

## Examples
To run the example scripts, you will need to have the spotless toolbox and MOSEK installed and added to your path. This version contains two examples: minimizing a quartic over the sphere and computing a set over which a network of Duffing oscillators is stable. The first problem solves
```
maximize        g
subject to      p(x) - g*(x'*x)^2 is SOS
```
to obtain a lower bound on `p(x)`. Because AnySOS produces feasible iterates, every iteration gives a lower bound on `p(x)` converging to the solution of the SDP. We compare this approach to two other solvers that trade speed for accuracy in the solution (DSOS and SDSOS from spotless) as well as solving the full SDP using MOSEK. These solvers only produce a feasible point in the limit of convergence, so they cannot quickly give a valid lower bound on `p(x)`.
