% Number of subintervals that we are going to create
N = 3000;
% Upper limit of the interval ([0,L])
L = 1;
% Step of the net
h = L / N;

% Net points
t_i = 0:h:L;

% Newton method constants
epsilon = 1 / 100;
num_iter = 50;

% Maximum lambda that we are going to evaluate
MAX_LAMBDA = 100;
lambda_h_slow = 1/100;
lambda_h_fast = 3/20;
low_limit_interval = 5;
upper_limit_interval = 70;

% Number of solutions that we are going to show in the same figure (the
% more, the less visible
NUMBER_SOLUTIONS = 15;

% Anulation interval
minimum_anulation = L / 2 - 1/10;
maximum_anulation = L / 2 + 1/10;
% Function a, the parameter
a = @(t)1 .* ((t < minimum_anulation) | (t > maximum_anulation));
degenerate = true;

% Discretize a (over all the net points)
a_eval = arrayfun(a,t_i);

% Different available colours for plotting solutions and bifurcation diag
colours = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];