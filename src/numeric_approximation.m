% Load all the constants
constants;

[u, bifurcation_values] = compute_solutions(h, N, L, t_i, a, a_eval, MAX_LAMBDA, lambda_h, epsilon, num_iter);

% Generate the bifurcation diagram
create_bifurcation_diagram(bifurcation_values, colours, 1);

% Generate the plot for the solutions
display_solutions(t_i, u, colours, NUMBER_SOLUTIONS, 2);

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

function [solutions_to_plot] = get_functions_to_plot(solutions, number_solutions)
% get_functions_to_plot function gets equally spaced functions from a range
% of solutions
% INPUT: 
%   - solutions: complete range of solutions
%   - number_solutions: maximum number of solutions that we want to select
% OUTPUT:
%   - solutions_to_plot: selected solutions
    step = (size(solutions, 1) - 1) / number_solutions;
    indexes = floor(1:step:size(solutions, 1) - 1);
    solutions_to_plot = zeros(number_solutions, size(solutions, 2));
    for i=1:number_solutions
        solutions_to_plot(i, :) = solutions(indexes(i),:);
    end
end

function display_solutions(t_i, solutions, colours, number_solutions, index)
% display_solutions function plots all the desired solutions in the same
% figure
% INPUT:
%   - t_i: points of the net
%   - solutions: all the solutions that the algorithm has obtained
%   - colours: all the available colours to select for the plots
%   - number_solutions: maximum number of solutions per eigenvalue
%   - index: number of the figure we want to plot
    number_rows = ceil(size(solutions, 2) / 2);
    figure(index);
    for i = 1:number_rows
        for j = 1:2
            % When there aren't more solutions to plot
            if 2 * (i - 1) + j > size(solutions, 2)
                break;
            end
            solutions_to_plot = get_functions_to_plot(solutions{1, 2 * (i - 1) + j}, number_solutions);
            % Create the subplot
            for k = 1:size(solutions_to_plot, 1)
                subplot(number_rows, 2, 2 * (i - 1) + j);
                plot(t_i, solutions_to_plot(k,:), colours(mod(k, size(colours, 2)) + 1))
                if 2 * (i - 1) + j == 1 % Positive solution
                    title('Soluciones positivas');
                elseif 2 * (i - 1) + j == 2 % 1 node solution (exclude s - plural)
                    title(sprintf('Soluciones con %d nodo', 2 * (i - 1) + j - 1));
                else
                    title(sprintf('Soluciones con %d nodos', 2 * (i - 1) + j - 1));
                end
                hold on;
            end
            hold off;
        end    
    end
end

function create_bifurcation_diagram(bifurcation_values, colours, index)
% create_bifurcation_diagram plots the global bifurcation diagram for the
% given solutions and lambdas
% INPUT:
%   - bifurcation_values: cell that contains for every row (eigenvalue), 
%           the lambdas in which we have obtained the solutions and their 
%           maximum values
%   - colours: all the available colours to select for the plots
%   - index: number of the figure we want to plot
    figure(index);
    for i = 1:size(bifurcation_values, 1)
        colour = colours(mod(i, size(colours, 2)) + 1);
        plot(bifurcation_values{i, 1}, bifurcation_values{i, 2}, colour);
        hold on;
        % Simmetric branch
        plot(bifurcation_values{i,1}, (-1).* bifurcation_values{i, 2}, colour);
        hold on;
    end
    
    % Configuration of the plot
    t = title("Diagrama de bifurcación global");
    % Locate both axis at the origin
    ax = gca;
    ax.XAxisLocation = 'origin';
    % Add x label and position it at the x-axis
    xlh = xlabel('\lambda');
    xlh.Position(2) = 0;
    
    % Add y label and position it at the top, without 90º rotation
    ylh = ylabel('||u||_\infty');
    ylh.Position(2) = t.Position(2)*1.02;
    ylh.Rotation = 0;
    
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

function [u, bifurcation_values] = compute_solutions(h, N, L, t_i, a, a_eval, limit_lambda, lambda_h, epsilon, num_iter)
% compute_solutions function obtains all the solutions for the problem that
% we are studying
% INPUT:
%   - h: step of the net
%   - N: number of subintervals
%   - L: upper limit of the interval
%   - t_i: points of the net
%   - a: problem parameter
%   - a_eval: discrete evaluation of a in t_i
%   - limit_lambda: maximum lambda for which we are going to obtain the
%           solutions
%   - lambda_h: step of lambda
%   - epsilon: maximum admisible error for the Newton method
%   - num_iter: maximum number of iterations for the Newton method before
%           considering that the method has diverged
% OUTPUT:
%   - u: computed solutions of diferent number of nodes. It is a cell-array
%           that contains, for each eigenvalue of the problem, a matrix, 
%           inside which, for every row, we have the values of the function 
%           evaluated in t_i
%   - bifurcation_values: maximum of each u, along the lambda for which
%           that function was obtained. It is a cell array that contains,
%           for each row, on the first position the values of lambda and on
%           the second position the maximum of each solution

    % 1 / h^2
    inverse_h_square = 1 / (h * h);

    % Obtain M (-u'' approximate)
    M = generate_M(N + 1, inverse_h_square);

    % Differential equation is transformed to algebraic equation
    F = @(lambda, x) -lambda * x + diag(a_eval) * x.^3 + M * x;
    
    % First eigenvalues that are limited by the MAX_LAMBDA constant
    eig = generate_eigenvalues_limited(inverse_h_square, N, limit_lambda);
    
    % Generate solutions for bifurcation branches from the bifurcation point in
    % our range (0, MAX_LAMBDA)
    for i = 1:size(eig, 2)
        % Initialize this iteration using crandall-rabinowitz theorem
        [u_0, lambda] = initialize_crandall_rabinowitz(i, L, h, t_i, a, eig(i), lambda_h);

        % Net for lambda
        lambdas = [eig(i) lambda:lambda_h:limit_lambda];
        % Zero solution is the first
        u_j = zeros(size(lambdas,2), N + 1);
        for j = 2:size(lambdas, 2)
            % The function we want to obtain the zeroes from is F but with a 
            % fixed lambda
            f = @(x) F(lambdas(j), x);
            jacobian_u_0 = M - lambdas(j) * eye(N + 1) + 3 * diag((diag(a_eval) * u_0.^2));
            [nu, ~, ~] = newton_modified_method(f, jacobian_u_0, u_0, epsilon, num_iter);
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
    end
end