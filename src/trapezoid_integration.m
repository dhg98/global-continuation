function [result] = trapezoid_integration(f, t_i, h)
% trapezoidIntegration calculates the approximate defined integral of a
% given function f, using the trapezoid method with subintervals
% INPUT:
%   - f: anonymous function for which we want to obtain the integral
%   - t_i: partition of the interval in which we want to calculate the
%   integral (given by i_limit = t_i[1], u_limit = t_i[size(t_i,2)])
%   - h: step of the partition
% OUTPUT:
%   - result: approximate value of the integral

    subst = f(t_i);
    result = (2 * sum(subst) - subst(1) - subst(t_i(size(t_i, 2)))) * h / 2;
end
