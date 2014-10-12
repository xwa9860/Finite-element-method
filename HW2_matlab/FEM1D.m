%{
This program 
- solve the 1D equation -d/dx[A_1(x)*du/dx] = f(x)
- with two boundary conditions, Dirichlet or Neumann
- using Gaussian quadrature for integration, order q
- uniform linear/quadratic/cubic elements, order p

***************************************************************************
INPUTS
***************************************************************************
N:   number of nodes. Integer  
A1Func:   fct A1(x) on the left side of the equation
fFunc:    f(x) on the right side of the equation
x0, xl:   Minimum and maximum coordinates of the problem domain.
v0, vl:   BC value at x0, BC value at xl
  

***************************************************************************
VARIABLES
***************************************************************************

Ne          : number of elements, Ne= N-1. Integer
coord       : array of coordinates of nodes. Array [N*1]. Can be changed to
                accomodate different sized elements
Q: number of GQ points used for integrals
GQpoint     : Coordinates of Gauss Quadrature points. Array [Q*1].
GQweight    : Weights of Gauss Quadrature points. Array [Q*1].
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
% Define the problem
% ************************************************************************
L = 1;
A1Func  =  @(x) -0.2;
fFunc  =  @(x) 256*sin(3/4*pi*x)*cos(16*pi*x);

xMin   =  0;
xMax   =  L;

bc_l_type = 1; %Dirichelet BC
bc_r_type = 2; %Neumann BC

u0    = 0;    % boundary value at xMin
tL    = -1;    % bourdary value at xMax

% ************************************************************************
% Define the mesh
% ************************************************************************

Ne =   10;
P = 1;
if P==1
    N=Ne+1;
%elseif P==2
end
  
coord = linspace(xMin, xMax, N);
matA = make_K_sparse_fast(1, Ne, coord, A1Func, fFunc, bc_l_type, u0, bc_r_type, tL);

%find the optimal number of element to achieve the error criteria
errorTol = 0.01;

%calculate the error of approximation  

error = calcError( Ne, A1Func, matA,coord);

while error > errorTol  && Ne < 1000
    Ne = 2* Ne;
    N=Ne+1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_sparse_fast(1, Ne, coord, A1Func, fFunc, bc_l_type, u0, bc_r_type, tL);
    error = calcError( Ne, A1Func, matA,coord);
end
display(Ne);
Ne_max = Ne;
Ne_min = Ne/2;

while Ne_max-Ne_min > 2
    Ne_Opt = ceil((Ne_max+Ne_min)/2);
    N=Ne_Opt+1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_sparse_fast(1, Ne_Opt, coord, A1Func, fFunc, bc_l_type, u0, bc_r_type, tL);
    error = calcError( Ne_Opt, A1Func, matA,coord);
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
[realU(e), du] = analyticalSolution(coord(e));
end

plot (coord, realU, 'r--');

grid on;
xlabel('x');
ylabel('u(x), u_N(x)');
title('Solution U(x)');
    
hold(gca, 'off');

% errorN = zeros(Ne_Opt+1, 1);
% 
% for e=1:Ne_Opt+1
% errorN(e) = analyticalSolution(coord(e), k)-matA(e);
% end
% 
% plot (coord, errorN, 'r--');


