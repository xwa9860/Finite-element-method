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

Ne =   1;  
P =2;

    
N = P*Ne +1;

coord = linspace(xMin, xMax, N);


matA = make_K_R(P, Ne, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
 % use backslash operator of MATLAB instead of inv()


% %find the optimal number of element to achieve the error criteria
errorTol = 0.010;
error =1;
%calculate the error of approximation  

while error > errorTol  && Ne<10000
    Ne = 2* Ne;
    N = P*Ne +1;
    coord = linspace(xMin, xMax, N);
    matA = make_K_R(P, Ne, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
    error = calcError(P, Ne, A1Func, dux, matA,coord);
    %display(error);
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
N_Opt = Ne_Opt*P +1;
display(N_Opt);
%display(Ne);
display(error);

% ************************************************************************
% Output the results and plot the solution
% ************************************************************************

% Defaults for this blog post
width = 5;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 30;      % Fontsize
lw = 5;      % LineWidth
msz = 9;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

% Plot the numerical solution, N_ksi =10 points per element
N_ksi =10;
index =1;
x_ksi_array = zeros(N_ksi*Ne_Opt,1);
u_ksi_numerical_array = zeros(N_ksi*Ne_Opt,1);
u_real_array = zeros(N_ksi*Ne_Opt,1);


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
        u_real_array(index) = ux(x_ksi);
        error_array(index) = abs(ux(x_ksi) -u);
        error_r_array(index) = abs((ux(x_ksi) -u)/u);
        index = index +1;
    end
end%for e = 1:NE_Opt

sol_fig = plot(x_ksi_array, u_ksi_numerical_array, 'o', x_ksi_array, u_real_array, 'r');
% 
% %realU = ux(coord);
% %plot (coord, realU);
% 
% set(gca,'FontSize',fsz);
% xlabel('x','FontSize', fsz);
% ylabel('u(x), u_N(x)','FontSize', fsz);
% %title('Solution for differential equation', 'FontSize', fsz);
% leg= legend('numerical', 'analytical', 'Location','southeast');
% set(leg,'FontSize',fsz);
% 
% set(gca,'XTick',-1:1); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%claculate error for different number of meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% P=1;
% Ne_Opt=1400;
% Ne_array = linspace(1, Ne_Opt, Ne_Opt);
% DOF_array = zeros(Ne_Opt, 1);
% h_array = zeros(Ne_Opt, 1);
% error_array = zeros(Ne_Opt, 1);
% for k = 1 :30: Ne_Opt
%     N = P * Ne_array(k) +1;
%     coord = linspace(xMin, xMax, N);
%     matA = make_K_R(P, k, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
%     error_array(k) = calcError(P, k, A1Func,dux, matA,coord);
%     DOF_array(k) = N;
%     h_array(k) = L/N;
% end
% 
% figure;
% error_fig = plot (h_array, error_array);
% 
% 
% 
% 
% hold all on;
% 
% P=2;
% Ne_Opt=100;
% Ne_array = linspace(1, Ne_Opt, Ne_Opt);
% DOF_array = zeros(Ne_Opt, 1);
% h_array = zeros(Ne_Opt, 1);
% error_array = zeros(Ne_Opt, 1);
% for k = 1 :2: Ne_Opt
%     N = P * Ne_array(k) +1;
%     coord = linspace(xMin, xMax, N);
%     matA = make_K_R(P, k, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
%     error_array(k) = calcError(P, k, A1Func,dux, matA,coord);
%     DOF_array(k) = N;
%     h_array(k) = L/N;
% end
% 
% error_fig = plot (h_array, error_array);
% 
% P=3;
% Ne_Opt=40;
% Ne_array = linspace(1, Ne_Opt, Ne_Opt);
% DOF_array = zeros(Ne_Opt, 1);
% h_array = zeros(Ne_Opt, 1);
% error_array = zeros(Ne_Opt, 1);
% for k = 1 :1: Ne_Opt
%     N = P * Ne_array(k) +1;
%     coord = linspace(xMin, xMax, N);
%     matA = make_K_R(P, k, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
%     error_array(k) = calcError(P, k, A1Func,dux, matA,coord);
%     DOF_array(k) = N;
%     h_array(k) = L/N;
% end
% 
% %error_fig = plot (DOF_array, error_array);
% error_fig = plot (h_array, error_array);
% set(gca,'FontSize',fsz);
% xlabel('element size','FontSize', fsz);
% ylabel('error','FontSize', fsz);
% %title('Evoluation of numerical error with number of elements', 'FontSize', fsz);
% leg= legend('P=1', 'P=2', 'P=3');
% %set(leg,'FontSize',fsz)
% 
% 
% % figure;
% % error_fig = plot (h_array, error_array);
% % set(gca,'FontSize',fsz);
% % xlabel('element sizes','FontSize', fsz);
% % ylabel('error','FontSize', fsz);
% % 
% 
