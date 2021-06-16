% Number of subintervals that we are going to create, less 1 (N + 2 points)
N = 2999;
% Upper limit of the interval ([0,L])
L = 1;
% Step of the net
h = L / (N + 1);

% Net points
t_i = 0:h:L;

% Newton method constants
big_epsilon = 2 / 10;
small_epsilon = 1 / 100;
num_iter = 100;

% Maximum lambda that we are going to evaluate
MAX_LAMBDA = 100;
lambda_h = [1/1000000, 1/1000, 3/20, 1/100];
limit_interval = [13/1000, 3, 70];

% Number of solutions that we are going to show in the same figure (the
% more, the less visible
NUMBER_SOLUTIONS = 15;

% Anulation interval
minimum_anulation = L / 2 - 7/20;
maximum_anulation = L / 2 + 7/20;
% Function a, the parameter
a = @(t)1 .* ((t < minimum_anulation) | (t > maximum_anulation));
degenerate = true;

% Discretize a (over all the net points)
a_eval = arrayfun(a,t_i);
a_eval(1) = [];
a_eval(end) = [];

% Different available colours for plotting solutions and bifurcation diag
colours = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];