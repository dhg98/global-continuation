function [x,fx,iters] = newton_modified_method(f, jacobian_x_0, x_0, epsilon, numIter)
% Newton Method computes the vectorial modified Newton algorithm to obtain 
% the root of a given function
%
% INPUT:
%   - f: anonymous function from which we want to obtain the root
%   - jacobian_x_0: Jf(x_0) (matrix)
%   - x_0: initial aproximation of the root (column)
%   - epsilon: method will halt when ||x_i - x_{i-1}||_2 < epsilon
%   ||x_{i-1}||_2
%           and ||f(x_i)||_2 < epsilon
%   - numIter: maximum number of iterations to consider the method has
%           diverged
% OUTPUT:
%   - x: last iteration point obtained by the algorithm (column)
%   - fx: evaluation of f(x)
%   - iters: number of iterations that the method has needed to get to a
%           good approximation
    
    % Initialization
    x0 = x_0;
    fx = f(x0);
    for i = 1:numIter
        % Next point
        x1 = x0 - jacobian_x_0\fx;
        
        % Evaluate function at the new point
        fx = f(x1);
        
        % End condition
        if norm(x1 - x0) < epsilon * norm(x0) && norm(fx) < epsilon
            break;
        end
        
        % Update variable for next iteration
        x0 = x1;
    end
    iters = i;
    x = x1;
end

