function [ux, dux] = analyticalSolution( x)
% ************************************************************************
% Calculate the analytical solution u(x) for a given x
% kcons is a constant in the equation
% ************************************************************************
ux = 5*((4087 * pi *x + 1536*sqrt(2) *x + ...
    137216*sin(61*pi*x/4)/61/pi -124928*sin(67*pi*x/4)/67/pi)...
    /4087/pi);
dux = 5*(1+512*(67*cos(61*pi*x/4) -...
    61*cos(67*pi*x/4))/4087/pi  ...
    -(512*(67*cos(61*pi/4) -...
    61*cos(67*pi/4))/4087/pi));
end

