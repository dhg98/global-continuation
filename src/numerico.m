% Load all the constants
constants;

% Obtain M (-u'' approximate)
M = generate_M(N + 1, inverse_h_square);

% Differential equation is transformed to algebraic equation
F = @(lambda, x) -lambda * x + diag(a_eval) * x.^3 + M * x;

% First eigenvalues that are limited by the MAX_LAMBDA constant
eig = generate_eigenvalues_limited(inverse_h_square, N, MAX_LAMBDA);
%bifurcation_points = 1;%size(eig, 2);

% u and bifurcation_values are uninitialized, so we need to clear them 
% before starting the loop
clear u;
clear bifurcation_values;

% Generate solutions for bifurcation branches from the bifurcation point in
% our range (0, MAX_LAMBDA)
for i = 1:size(eig, 2)
    % Initialize this iteration using crandall-rabinowitz theorem
    [u_0, lambda] = initialize_crandall_rabinowitz(i, L, h, t_i, a, eig(i), lambda_h);
    
    % Net for lambda
    lambdas = [eig(i) lambda:lambda_h:MAX_LAMBDA];
    % Zero solution is the first
    u_j = zeros(size(lambdas,2), N + 1);
    for j = 2:size(lambdas, 2)
        % The function we want to obtain the zeroes from is F but with a 
        % fixed lambda
        f = @(x) F(lambdas(j), x);
        jacobian_u_0 = M - lambdas(j) * eye(N + 1) + 3 * diag((diag(a_eval) * u_0.^2));
        [nu, ~, iters] = newton_modified_method(f, jacobian_u_0, u_0, epsilon, num_iter);
        % Next u_0 is from where we finished in our last iteration
        u_0 = nu;
        u_j(j,:) = nu;
        %plot(t_i, nu, colours(mod(j, size(colours, 2)) + 1));
        %hold on;
    end
    %hold off;
    u{i} = u_j;
    bifurcation_values{i,1} = lambdas;
    bifurcation_values{i,2} = max(u{i},[],2);
    % create_bifurcation_branch(squeeze(u(i,:,:)), lambdas, 'r');
end
create_bifurcation_diagram(bifurcation_values, colours);

function [M] = generate_M(n, inverse_h_square)
% generateM generates the finite differences matrix of -u''
% INPUT:
%   - n: size of the matrix
%   - inverse_h_square: value of 1/h^2
% OUTPUT:
%   - M: generated finite differences matrix
    M = inverse_h_square * 2 * eye(n);
    for i = 1:n-1
        M(i + 1, i) = -inverse_h_square;
        M(i, i + 1) = -inverse_h_square;
    end
end

function [eig] = generate_ith_eigenvalue(i, inverse_h_square, sizeM)
% generateEigenvalues calculates the eigenvalue i of the matrix M
% defined previously, provided that M is a Toeplitz matrix, so that the 
% eigenvalue follows a fixed formula
% INPUT:
%   - i: index of eigenvalue we want to calculate (between 1 and sizeM both
%   included
%   - inverse_h_square: value of 1/h^2
%   - sizeM: size of the matrix M
% OUTPUT:
%   - eig: list of eigenvalues that have been obtained

    eig = 2 * inverse_h_square * (1 - cos((i * pi)/ (sizeM + 1)));
end

function [eigs] = generate_eigenvalues_limited(inverse_h_square, N, limit_lambda)
% generate_eigenvalues_limited obtains the first eigenvalues based on a 
% limit on lambda (that is, we get all the eigenvalues that are under that
% limit. In order to do so, we will use that M is Toeplitz
% INPUT:
%   - inverse_h_square: inverse of the square of the net step
%   - N: number of subintervals
%   - limit_lambda: maximum lambda that we want to explore
% OUTPUT:
%   - eigs: approximate eigenvalues that are less than limit_lambda

    i = 1;
    while 1
        eigenvalue = generate_ith_eigenvalue(i, inverse_h_square, N + 1);
        if eigenvalue >= limit_lambda || i == N + 2
            break;
        end
        eigs(i) = eigenvalue;
        i = i + 1;
    end
end

function [s_phi_0] = generate_s_phi_0(n, L, t_i, s)
% generate_s_phi_0 obtains the first discrete function for the newton
% method, which is s * phi_0, where phi_0(t) = sin(n * pi * t / L)
% INPUT:
%   - n: the parameter for which we want to generate the autofunction
%   - L: the upper limit of the interval
%   - t_i: points of the net
%   - s: value of the parameter
% OUTPUT:
%   - s_phi_0: calculation of s * phi_0(t) in all the grid points (column)
    
    f = @(t) sin(n .* pi .* t ./ L);
    s_phi_0 = (s .* f(t_i)).';
end

function create_bifurcation_diagram(bifurcation_values, colours)
% create_bifurcation_diagram plots the global bifurcation diagram for the
% given solutions and lambdas
% INPUT:
%   - bifurcation_values: cell that contains for every row (eigenvalue), 
%           the lambdas in which we have obtained the solutions and their 
%           maximum values
    for i = 1:size(bifurcation_values, 1)
        plot(bifurcation_values{i, 1}, bifurcation_values{i, 2}, colours(mod(i, size(colours, 2)) + 1));
        hold on;
    end
    title("Diagrama de bifurcación global");
    hold off;
end

function [u_0, lambda] = initialize_crandall_rabinowitz(i, L, h, t_i, a, eig, lambda_h)
% initialize_crandall function initializes the approximate method for
% the i-th eigenvalue (branch that will start at this eigenvalue) using the
% Crandall-Rabinwitzt theorem
% INPUT:
%   - i: iteration
%   - L: upper limit
%   - h: step of the net
%   - t_i: points of the net
%   - a: parameter of the problem
%   - eig: eigenvalue that we are taking into account
%   - lambda_h: step for the lambda (distance to 'exact' eigenvalue)
% OUTPUT:
%   - u_0: s * phi_0 whose existence is assured by Crandall-Rabinowitz 
%           theorem
%   - lambda: lambda_0 + s^2 lambda_2 that is obtained based on the
%           Crandall-Rabinowitz theorem

    % Calculation of lambda_2
    phi_a = @(t) a(t) .* (sin(i .* pi .* t ./ L)).^4;
    lambda_2 = trapezoid_integration(phi_a, t_i, h) * 2 / L;
    % Value of the parameter that we get from the Crandall-Rabinowitz theorem
    s = sqrt(((i * pi / L)^2 + lambda_h - eig) / lambda_2);
    % Generate initial lambda
    lambda = eig + s^2 * lambda_2;
    % Generate initial u_0 (s * phi_0) (column)
    u_0 = generate_s_phi_0(i, L, t_i, s);
end