# Global continuation of nodal solutions

This project is part of a final degree project made by Daniel Herranz at University Complutense of Madrid for the Double Bachelor in Computer Science and Mathematics.

The intention of this project is to represent nummerically all the solutions of the boundary value problem `-u'' = λu - a(t)u^3` in `[0,L]`, with `u(0) = 0`, `u(L) = 0`, where `L > 0` and `a(t)` is a non-negative function which is positive in some point, and the global bifurcation diagram of these solutions. However, it is possible to adapt the functioning of the program to a different differential equation by using the same finite difference method, but applied to different functions.

This project is the culmination of a exhaustive study of the properties the solutions must have, that begins studying the Netwon method and the Implicit function theorem and ends with the Crandall-Rabinowitz theorem, which ensures a bifurcation branch from points (aka eigenvalues) of the form `(nπ/L)^2`, for all `n≥1`, and which is used to initialize the branch for every eigenvalue.

## Installation
In order to make the project work, you need to have MATLAB installed, and a license for it. After that step is done, simply run the `numeric_approximation.m` file and you will be ready to go. Only modify the `constants.m` file if you know what you are doing, because it may cause errors in the execution of the approximation program, or slow down the processing (it can even cause the Newton method to diverge).

## Output
As it has already been said, this program tries to plot all the different solutions with a limited number of nodes (that is, interior zeroes) for a limited λ, and also computes the global bifurcation diagram (plots λ vs. the infinity norm of the computed functions, that is `||u||_\infty`) based on the solutions that have been obtained.

We have computed the profile of the solutions for two different cases, one <b>classical</b>, that can be analyzed theorically without the use of Crandall-Rabinowitz theorem, by only using phase diagram theory, and the other one <b>degenerate</b>, that uses more complex mathematics.

### Classical
In this case, we are obtaining the results for the problem `-u'' = λu - u^3` in `[0,1]`, with `u(0) = 0`, `u(1) = 0`. Note that here, `a` is constant, and equal to 1, which satisfies the condition we have stated previously (`a` is non-negative and non 0 constant). In the following images, you can see the profile of the different solutions, together with the bifurcation diagram.

![classical-grid-solutions](https://github.com/dhg98/global-continuation/blob/main/res/classic/grid_solutions.png)

![classical-global-bifurcation-diagram](https://github.com/dhg98/global-continuation/blob/main/res/classic/global_bifurcation_diagram.png)

### Degenerate crossed bifurcation diagram
In this case, we are obtaining the results for the problem `-u'' = λu - a(t)u^3` in `[0,1]`, with `u(0) = 0`, `u(1) = 0`, where `a` is the caracteristic function over the set `[0, 0.4) \cup (0.6, 1]`. Note that in this case we can't use phase diagram technique, so we need to strongly use Crandall-Rabinowitz theorem. In the following image you can see the profile of the different solutions, together with the bifurcation diagram.

![degenerate1-grid-solutions](https://github.com/dhg98/global-continuation/blob/main/res/degenerate_crossed/grid_solutions.jpg)

![degenerate1-global-bifurcation-diagram](https://github.com/dhg98/global-continuation/blob/main/res/degenerate_crossed/global_bifurcation_diagram.jpg)

### Degenerate uncrossed bifurcation diagram
In this case, we are obtaining the results for the problem `-u'' = λu - a(t)u^3` in `[0,1]`, with `u(0) = 0`, `u(1) = 0`, where `a` is the caracteristic function over the set `[0, 0.15) \cup (0.85, 1]`. Note that in this case we can't use phase diagram technique, so we need to strongly use Crandall-Rabinowitz theorem. In the following image you can see the profile of the different solutions, together with the bifurcation diagram.

![degenerate2-grid-solutions](https://github.com/dhg98/global-continuation/blob/main/res/degenerate_not_crossed/grid_solutions.png)

![degenerate2-global-bifurcation-diagram](https://github.com/dhg98/global-continuation/blob/main/res/degenerate_not_crossed/global_bifurcation_diagram.png)