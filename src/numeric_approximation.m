% Load all the constants
constants;

% 1 / h^2
inverse_h_square = 1 / (h * h);
% Obtain all the eigenvalues that we are going to compute
eig = generate_eigenvalues_limited(inverse_h_square, N, MAX_LAMBDA);
if exist('u','var') == 0
    u = {};
    bifurcation_values = {};
    for i = 1:size(eig, 2)
        u{i} = [];
        bifurcation_values{i, 1} = [];
        bifurcation_values{i, 2} = [];
    end
end

err_indexes = [];
for i = 1:size(eig,2)
    if i > size(u, 2)
        u{i} = [];
        bifurcation_values{i, 1} = [];
        bifurcation_values{i, 2} = [];
    end
    if size(u{i}, 1) == 0
        % Eigen value i doesn't have its solutions loaded
        [u_aux, bif_aux, err] = load_solutions_eigenvalue(degenerate, i);
        u{i} = u_aux;
        bifurcation_values{i, 1} = bif_aux{1};
        bifurcation_values{i, 2} = bif_aux{2};
        clear u_aux  bif_aux;
        if err
            err_indexes = [err_indexes, i];
        end
    end
end

% There are files with errors. We need to compute the solutions
if size(err_indexes, 2) > 0
    % Obtain M (-u'' approximate)
    M = generate_M(N + 1, inverse_h_square);

    % Differential equation is transformed to algebraic equation
    F = @(lambda, x) -lambda * x + diag(a_eval) * x.^3 + M * x;
    
    for i = 1:size(err_indexes, 2)
        [u_aux, bif_aux] = compute_branch(h, N, L, t_i, a, a_eval, i, eig(i), lambda_h_slow, lambda_h_fast, low_limit_interval, upper_limit_interval, epsilon, num_iter, maximum_anulation, minimum_anulation, degenerate, err_indexes);
        u{i} = u_aux;
        bifurcation_values{i, 1} = bif_aux{1};
        bifurcation_values{i, 2} = bif_aux{2};
        clear u_aux  bif_aux;
    end
    clear F M inverse_h_square eig;
end

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
        leg{i} = sprintf("Bifurcación desde autovalor %0.5f", bifurcation_values{i, 1}(1));
        hold on;
    end
    legend(leg);
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

function [u, bifurcation_values] = compute_branch(h, N, L, t_i, a, a_eval, i, eig, lambda_h_slow, lambda_h_fast, low_limit_interval, upper_limit_interval, epsilon, num_iter, maximum_anulation, minimum_anulation, degenerate)
% compute_branch function obtains all the solutions for the problem that
% we are studying for a specific eigenvalue branch
% INPUT:
%   - h: step of the net
%   - N: number of subintervals
%   - L: upper limit of the interval
%   - t_i: points of the net
%   - a: problem parameter
%   - a_eval: discrete evaluation of a in t_i
%   - i: iteration for which we are obtaining the branch
%   - eig: eigenvalue for which we are obtaining the branch
%   - lambda_h_slow: slow step of lambda (for the first period and for the
%       last one
%   - lambda_h_fast: fast step of lambda (for the middle period)
%   - low_limit_interval: distance from the eigenvalue to use slow step
%   - upper_limit_interval: distance from the maximum possible eigenvalue
%       in degenerate cases
%   - minimum_anulation: low limit interval of anulation of a (for
%       degenerate cases)
%   - maximum_anulation: upper limit interval of anulation of a (for
%       degenerate cases)
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
    
    lambda_step = lambda_h_slow;
    % Initialize this iteration using crandall-rabinowitz theorem
    [u_0, lambda] = initialize_crandall_rabinowitz(i, L, h, t_i, a, eig, lambda_step);
    % Net for lambda
    lambdas = [eig, lambda];
    j = 2;
    % Zero solution is the first
    u_j = zeros(1, N + 1);
    max_lambda_value = (i * pi / (maximum_anulation-minimum_anulation))^2;
    while 1
        % The function we want to obtain the zeroes from is F but with a 
        % fixed lambda
        f = @(x) F(lambdas(j), x);
        jacobian_u_0 = M - lambdas(j) * eye(N + 1) + 3 * diag((diag(a_eval) * u_0.^2));
        [nu, ~, iters] = newton_modified_method(f, jacobian_u_0, u_0, epsilon, num_iter);
        disp(iters);
        if  iters == num_iter % Max iterations -> method has diverged
            lambdas(j) = [];
            break;
        end

        if lambdas(j) - lambdas(1) > low_limit_interval
            lambda_step = lambda_h_fast;
        end
        if max_lambda_value - lambdas(j) < upper_limit_interval
            lambda_step = lambda_h_slow;
        end

        % Next u_0 is from where we finished in our last iteration
        u_0 = nu;
        u_j(j,:) = nu;

        % Generate new lambda based on iterations in last newton execution
        new_lambda = lambdas(j) + lambda_step / iters;
        if new_lambda > max_lambda_value
            break;
        else
            lambdas(j + 1) = new_lambda;
        end
        j = j + 1;
    end
    u{i} = u_j;
    if degenerate
        extra = '_degenerate';
    else
        extra = '';
    end
    write_functions_in_file(sprintf('../data/data%s%d.txt', extra, i), u_j);
    bifurcation_values{i,1} = lambdas;
    bifurcation_values{i,2} = max(u{i},[],2);
    write_bifurcation_diag_in_file(sprintf('../data/bifurcation%s%d.txt', extra, i), lambdas, bifurcation_values{i,2});
end

function write_functions_in_file(filename, u)
% write_functions_in_file function saves all the functions computed at u
% variable in a file
% INPUT:
%   - filename: name that we are going to use to save the file. If it
%   exists, it will be overwritten
%   - u: (discrete) functions that will be saved in the file
    temp = u;
    save(filename, 'temp', '-ascii', '-double');
end

function write_bifurcation_diag_in_file(filename, lambdas, values)
% write_bifurcation_diag_in_file function saves all the functions computed
% at (lambdas, values) variables in a file
% INPUT:
%   - filename: name that we are going to use to save the file. If it
%   exists, it will be overwritten
%   - lambdas: array of lambda values that we have used to compute all the
%       different solutions
%   - values: array of values that match with the maximum of the function
%       computed for the lambda that has the same index
    if exist(filename, 'file')==2
        delete(filename);
    end
    for i = 1:size(lambdas, 2)
        temp = [lambdas(i), values(i)];
        save(filename, 'temp', '-append', '-ascii', '-double');
    end
end

function [u, bifurcation, err] = load_solutions_eigenvalue(degenerate, eigenvalue)
% load_solutions_eigenvalue function loads the solutions of a specific
% eigenvalue given by its index. It will try to find (depending on if we
% are considering the degenerate case or the classic case) the two files 
% that have the data for this case, and load its content to u and
% bifurcation variables. It will check for errors during opening the files,
% because that would mean there has been an error.
% INPUT:
%   - degenerate: boolean that indicates if we want to load a file that is
%   degenerated or not
%   - eigenvalue: index of the eigenvalue we are interested in
% OUTPUT:
%   - u: solutions loaded, or empty if there was an error
%   - bifurcation: bifurcation branch loaded, or empty if there was an
%       error
%   - err: boolean that indicates if there was an error or not
    ls = dir('../data');
    if (degenerate) 
        solutions_start = sprintf("data_degenerate%d", eigenvalue);
        eigenvalue_start = sprintf("bifurcation_degenerate%d", eigenvalue);
    else
        solutions_start = sprintf("data_classic%d", eigenvalue);
        eigenvalue_start = sprintf("bifurcation_classic%d", eigenvalue);
    end
    solutions_filename = "";
    bifurcation_filename = "";
    for i = 1:size(ls,1)
        if ~ls(i).isdir
            if startsWith(ls(i).name, solutions_start)
                solutions_filename = ls(i).name;
            end
            if startsWith(ls(i).name, eigenvalue_start)
                bifurcation_filename = ls(i).name;
            end
        end
    end
    
    % Solution
    bifurcation{1} = [];
    bifurcation{2} = [];
    u = [];
    solutions_filename = sprintf('../data/%s', solutions_filename);
    bifurcation_filename = sprintf('../data/%s', bifurcation_filename);
    fid1 = fopen(solutions_filename);
    fid2 = fopen(bifurcation_filename);
    err = 0;
    if fid1 == -1 || fid2 == -1
        err = 1;
    end
    if err == 0
        while true
            thisline1 = fgetl(fid1);
            if ~ischar(thisline1); break; end  % End of file
            % Both files have the same length
            thisline2 = fgetl(fid2);
            solution_values = strsplit(thisline1);
            bifurcation_values = strsplit(thisline2);
            % Erase first tab
            solution_values(1) = [];
            bifurcation_values(1) = [];
            u(size(u, 1) + 1,:) = str2double(solution_values);
            bifurcation{1}(size(bifurcation{1}, 2) + 1) = str2double(bifurcation_values(1));
            bifurcation{2}(size(bifurcation{2}, 1) + 1, 1) = str2double(bifurcation_values(2));
        end
        fclose('all');
    end
end