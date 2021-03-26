# Global continuation of nodal solutions

This project is part of a final degree project made by Daniel Herranz at University Complutense of Madrid for the Double Bachelor in Computer Science and Mathematics.

The intention of this project is to represent nummerically all the solutions of the boundary value problem `-u'' = λu - a(t)u^3` in `[0,L]`, with `u(0) = 0`, `u(L) = 0`, where `L > 0` and `a(t)` is a non-negative function which is positive in some point, and the global bifurcation diagram of these solutions. However, it is possible to adapt the functioning of the program to a different differential equation by using a similar strategy.

This project is the culmination of a exhaustive study of the properties the solutions must have, that begins studying the Netwon method and the Implicit function theorem and ends with the Crandall-Rabinowitz theorem, which ensures a bifurcation branch from points (aka eigenvalues) of the form `(nπ/L)^2`, for all `n≥1`, and which is used to initialize the branch for every eigenvalue. 

## Installation
In order to make the project work, you need to have MATLAB installed, and a license for it. After that step is done, simply run the `numeric_approximation.m` file and you will be ready to go. Only modify the `constants.m` file if you know what you are doing, because it may cause errors in the execution of the approximation program, or slow down the processing (it can even cause the Newton method to diverge).

## Output
As it has already been said, this program tries to plot all the different solutions with a limited number of nodes (that is, interior zeroes) for a limited λ, and also computes the global bifurcation diagram based on the solutions that have been obtained. Here you can see some examples of a few functions that the program is obtaining for the problem `-u'' = λu - u^3` in `[0,1]`, with `u(0) = 0`, `u(1) = 0`.

1. Positive solution
![positive-solution](https://github.com/dhg98/global-continuation/blob/main/res/positive_solution.jpg)

2. 1-node solution
![1-node-solution](https://github.com/dhg98/global-continuation/blob/main/res/1_node_solution.jpg)

3. 2-node solution
![2-node-solution](https://github.com/dhg98/global-continuation/blob/main/res/2_node_solution.jpg)

4. Global bifurcation diagram
![global-bifurcation-diagram](https://github.com/dhg98/global-continuation/blob/main/res/global_bifurcation_diagram.jpg)
