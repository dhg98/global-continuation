% Number of subintervals that we are going to create
N = 1000;
% Upper limit of the interval ([0,L])
L = 1;
% Step of the net
h = L / N;
h_square = h * h;
inverse_h_square = 1 / h_square;

% Net points
t_i = 0:h:L;

% Newton method constants
epsilon = 1 / 100;
num_iter = 20;

% Maximum lambda that we are going to evaluate
MAX_LAMBDA = 50;
lambda_h = 1/100;

% Function a, the parameter
a = @(t)1;

% Discretize a (over all the net points)
a_eval = arrayfun(a,t_i);

% Different available colours for plotting solutions and bifurcation diag
colours = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];