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

Ne  =   10;
N   =   Ne+1;    

k=2;
L=1;
A1Func  =  '-0.2';    % These can be functions of x, e.g. 'sin(x)+ x^2'
fFunc  =  'k^2*sin(pi*k*x/1)+ 2*x';%'k^2*sin(pi*k*x/L)+ 2*x';

xMin   =  0;
xMax   =  L;

u0     = 0;    % boundary value at xMin
uL     = 1;    % bourdary value at xMax


coord = linspace(xMin, xMax, N);

matA = make_K_sparse_fast(N,k, coord, A1Func, fFunc, u0, uL);
 % use backslash operator of MATLAB instead of inv()



%find the optimal number of element to achieve the error criteria
errorTol = 0.01;

%calculate the error of approximation  

error = calcError( Ne, A1Func, matA,k,coord);

while error > errorTol  
    Ne = 2* Ne;
    N=Ne+1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_sparse_fast(N,k, coord, A1Func, fFunc, u0, uL);
    error = calcError( Ne, A1Func, matA,k,coord);
end
display(Ne);
Ne_max = Ne;
Ne_min = Ne/2;

while Ne_max-Ne_min > 2
    Ne_Opt = ceil((Ne_max+Ne_min)/2);
    N=Ne_Opt+1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_sparse_fast(N,k, coord, A1Func, fFunc, u0, uL);
    error = calcError( Ne_Opt, A1Func, matA,k,coord);
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

% Plot the nodal values as circles.
plot(coord, matA, 'o', 'MarkerFaceColor', 'g');
hold(gca, 'on');

% Plot the solution over each element as lines. Straight lines are enough
% for linear elements.
for e = 1:Ne_Opt
   plot (coord, matA, 'g');  
end

realU = zeros(Ne_Opt, 1);

for e=1:Ne_Opt+1
[realU(e), du] = analyticalSolution(coord(e), k);
end

plot (coord, realU, 'r--');

grid on;
xlabel('x');
ylabel('u(x), u_N(x)');
title('Solution U(x)');
    
hold(gca, 'off');

errorN = zeros(Ne_Opt+1, 1);

for e=1:Ne_Opt+1
errorN(e) = analyticalSolution(coord(e), k)-matA(e);
end

plot (coord, errorN, 'r--');


