%{
This program 
- solve the 1D equation -d/dx[A_1(x)*du/dx] + A_2(x)u = f(x)
- with two Dirichlet boundary conditions.
- using 3 point Gaussian quadrature for integration
- uniform linear elements

***************************************************************************
INPUTS
***************************************************************************
N                    :   number of nodes. Integer  

A1Func, A2Func     :   fcts in the left side of the equation
fFunc                :   f(x) in the right side of the equation
xMin, xMax           :   Minimum and maximum coordinates of the problem domain.
u0, ul               :   BC u(0), BC u(L)
  

***************************************************************************
VARIABLES
***************************************************************************

Ne          : number of elements, Ne= N-1. Integer
coord       : array of coordinates of nodes. Array [N*1]. Can be changed to
                accomodate different sized elements
GQpoint     : Coordinates of Gauss Quadrature points. Array [3*1].
GQweight    : Weights of Gauss Quadrature points. Array [3*1].
S           : Linear shape functions evaluated at GQ points. 
dS          : Derivatives of shape functions evaluated at GQ points. 
matK        : stiffness matrix, sparse matrix N*N
matF        : Global force vector[1*N]. 
matA        : result matrix

%}

commandwindow
clc;          % Clear the command window
clear all;    % Clear the memory, delete previously stored data
close all;    % Close previously opened figure windows


% ************************************************************************
% Define variables
% ************************************************************************

A1Func  = @(x) -0.2;    % These can be functions of x, e.g. 'sin(x)+ x^2'
fFunc  = @(x) 256*sin(3/4*pi*x)*cos(16*pi*x);%'k^2*sin(pi*k*x/L)+ 2*x';

L=1;
xMin   =  0;
xMax   =  L;

type_bc_l=1;
v0     = 0;    % boundary value at xMin
type_bc_r=2;
vl     = -1;    % bourdary value at xMax

%analytical solution for the problem
ux =@(x) 5*((4087 * pi *x + 1536*sqrt(2) *x + ...
    137216*sin(61*pi*x/4)/61/pi -124928*sin(67*pi*x/4)/67/pi)...
    /4087/pi);

dux =@(x) 5*(512*(67*cos(61*pi*x/4) -...
    61*cos(67*pi*x/4))/4087/pi +1 ...
    -(512*(67*cos(61*pi/4) -...
    61*cos(67*pi/4))/4087/pi));

Ne =   10;  
P = 3;
N = P*Ne +1;

coord = linspace(xMin, xMax, N);


matA = make_K_R(P, Ne, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
 % use backslash operator of MATLAB instead of inv()


%find the optimal number of element to achieve the error criteria
errorTol = 0.01;

%calculate the error of approximation  

error = calcError(P, Ne, A1Func, dux, matA,coord);

while error > errorTol  && Ne<10000
    Ne = 2* Ne;
    N = P*Ne +1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_R(P, Ne, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
    error = calcError(P, Ne, A1Func, dux, matA,coord);
end
display(Ne);
Ne_max = Ne;
Ne_min = Ne/2;

while Ne_max-Ne_min > 2
    Ne_Opt = ceil((Ne_max+Ne_min)/2);
    N = P*Ne_Opt +1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_R(P, Ne_Opt, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
    error = calcError(P, Ne_Opt, A1Func,dux, matA,coord);
    if error > errorTol
        Ne_min = Ne_Opt;
        %display(Ne_min);
    else
        Ne_max = Ne_Opt;
        %display(Ne_max);
    end
end
display(Ne_Opt);
display(error);

% ************************************************************************
% Output the results and plot the solution
% ************************************************************************

% Plot the numerical solution, N_ksi =10 points per element
N_ksi =10;
index =1;
x_ksi_array = zeros(N_ksi*Ne_Opt,1);
u_ksi_numerical_array = zeros(N_ksi*Ne_Opt,1);

for e=1:Ne_Opt
    Nl= P*(e-1)+1;
    Nr= P*e+1;
    coordElement = coord(Nl:Nr);
    matA_elem=matA(Nl:Nr);
    ksi_array = linspace(-1, 1, N_ksi);
    for k = 1:1:N_ksi
        [du, u, x_ksi] = postprocessing (ksi_array(k), matA_elem, coordElement, P);
        x_ksi_array(index) = x_ksi;
        u_ksi_numerical_array(index) = u;
        index = index +1;
    end
end%for e = 1:NE_Opt

plot(x_ksi_array, u_ksi_numerical_array, 'o');
hold all;

realU = ux(coord);
% for e=1:Ne_Opt+1
% [realU(e), du] = analyticalSolution(coord(e));
% end

plot (coord, realU, '--');
%plot(coord, matA, 'o');
grid on;
xlabel('x');
ylabel('u(x), u_N(x)');
title('Solution for differential equation');
legend('numerical', 'analytical');
hold(gca, 'off');
% 
% errorN = zeros(Ne_Opt+1, 1);
% 
% for e=1:Ne_Opt+1
% errorN(e) = analyticalSolution(coord(e), k)-matA(e);
% end
% 
% plot (coord, errorN, 'r--');


