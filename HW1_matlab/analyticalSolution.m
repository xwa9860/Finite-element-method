function [ux, dux] = analyticalSolution( x, k )
% ************************************************************************
% Calculate the analytical solution u(x) for a given x
% kcons is a constant in the equation
% ************************************************************************
ux = -5/(pi)^2*sin(pi * k *x) + 5/3*x^3 - 2/3*x;

dux = -5*k/(pi)*cos(pi * k *x) + 5*x^2 - 2/3;
end

