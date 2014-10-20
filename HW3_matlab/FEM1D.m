%{
This program 
- solve the 1D equation d/dx[A_1(x)*du/dx] = f(x)
- using Gaussian quadrature for integration
- uniform linear elements
%}

commandwindow
clc;          % Clear the command window
clear all;    % Clear the memory, delete previously stored data
close all;    % Close previously opened figure windows

% ************************************************************************
% Define the equation
% ************************************************************************
fFunc  = @(x) 256*sin(3/4*pi*x)*cos(16*pi*x);%'k^2*sin(pi*k*x/L)+ 2*x';

L=1;
xMin   =  0;
xMax   =  L;

type_bc_l=1;
v0     = 0;    % boundary value at xMin
type_bc_r=2;
vl     = 1;    % bourdary value at xMax

%Ne_array= [100, 1000, 10000];
%for k = 1:3
Ne =   100;  
P =1;
N = P*Ne +1;
coord = linspace(xMin, xMax, N);
[KE_table, R]= make_K_R(P, Ne, coord, fFunc, type_bc_l, v0, type_bc_r, vl);
[matA, error] = CG_solver_KE_table(KE_table, R, 1e-8);
%matA = K\R;

%ux_arr_2 =zeros(N, 1);
ux_arr =zeros(N, 1);
for i=1:N
    ux_arr(i) = ana_sol_calculator(coord(i));
    %ux_arr_2(i) = ux_2(coord(i));
end
plot(coord, matA,'o');
hold all on;
optimized_plot(coord, ux_arr);
%optimized_plot(coord, ux_arr_2)

legend('numerical solution', 'analytical solution')

% % %find the optimal number of element to achieve the error criteria
% errorTol = 0.010;
% error =1;
% %calculate the error of approximation  
% 
% %while error > errorTol  && Ne<100
%     Ne = 2* Ne;
%     N = P*Ne +1;
%     coord = linspace(xMin, xMax, N);
%     [K, R]= make_K_R(P, Ne, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%     matA = CG_solver(K, R);
%     error = calcError(P, Ne, dux, matA,coord);
% %end
% display(Ne);
% Ne_max = Ne;
% Ne_min = Ne/2;
% display(error)
% 
% % while Ne_max-Ne_min > 2
% %     Ne_Opt = ceil((Ne_max+Ne_min)/2);
% %     N = P*Ne_Opt +1;
% %     coord = linspace(xMin, xMax, N);
% %     matA = make_K_R(P, Ne_Opt, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
% %     error = calcError(P, Ne_Opt, dux, matA,coord);
% %     if error > errorTol
% %         Ne_min = Ne_Opt;
% %         %display(Ne_min);
% %     else
% %         Ne_max = Ne_Opt;
% %         %display(Ne_max);
% %     end
% % end
% % display(Ne_Opt);
% % N_Opt = Ne_Opt*P +1;
% % display(N_Opt);
% % %display(Ne);
% % display(error);
% % 
% % % ************************************************************************
% % % Output the results and plot the solution
% % % ************************************************************************
% 
% % % Plot the numerical solution, N_ksi =10 points per element
% % N_ksi =10;
% % index =1;
% % x_ksi_array = zeros(N_ksi*Ne_Opt,1);
% % u_ksi_numerical_array = zeros(N_ksi*Ne_Opt,1);
% % u_real_array = zeros(N_ksi*Ne_Opt,1);
% 
% 
% % for e=1:Ne_Opt
% %     Nl= P*(e-1)+1;
% %     Nr= P*e+1;
% %     coordElement = coord(Nl:Nr);
% %     matA_elem=matA(Nl:Nr);
% %     ksi_array = linspace(-1, 1, N_ksi);
% %     for k = 1:1:N_ksi
% %         [du, u, x_ksi] = postprocessing (ksi_array(k), matA_elem, coordElement, P);
% %         x_ksi_array(index) = x_ksi;
% %         u_ksi_numerical_array(index) = u;
% %         u_real_array(index) = ux(x_ksi);
% %         error_array(index) = abs(ux(x_ksi) -u);
% %         error_r_array(index) = abs((ux(x_ksi) -u)/u);
% %         index = index +1;
% %     end
% % end%for e = 1:NE_Opt
% % 
% % sol_fig = plot(x_ksi_array, u_ksi_numerical_array, 'o', x_ksi_array, u_real_array, 'r');
% % 
% % %realU = ux(coord);
% % %plot (coord, realU);
% % 
% % set(gca,'FontSize',fsz);
% % xlabel('x','FontSize', fsz);
% % ylabel('u(x), u_N(x)','FontSize', fsz);
% % %title('Solution for differential equation', 'FontSize', fsz);
% % leg= legend('numerical', 'analytical', 'Location','southeast');
% % set(leg,'FontSize',fsz);
% % 
% % set(gca,'XTick',-1:1); %<- Still need to manually specific tick marks
% % set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %claculate error for different number of meshes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% % P=1;
% % Ne_Opt=1400;
% % Ne_array = linspace(1, Ne_Opt, Ne_Opt);
% % DOF_array = zeros(Ne_Opt, 1);
% % h_array = zeros(Ne_Opt, 1);
% % error_array = zeros(Ne_Opt, 1);
% % for k = 1 :30: Ne_Opt
% %     N = P * Ne_array(k) +1;
% %     coord = linspace(xMin, xMax, N);
% %     matA = make_K_R(P, k, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
% %     error_array(k) = calcError(P, k, A1Func,dux, matA,coord);
% %     DOF_array(k) = N;
% %     h_array(k) = L/N;
% % end
% % 
% % figure;
% % error_fig = plot (h_array, error_array);
% % 
% % 
% % 
% % 
% % hold all on;
% % 
% % P=2;
% % Ne_Opt=100;
% % Ne_array = linspace(1, Ne_Opt, Ne_Opt);
% % DOF_array = zeros(Ne_Opt, 1);
% % h_array = zeros(Ne_Opt, 1);
% % error_array = zeros(Ne_Opt, 1);
% % for k = 1 :2: Ne_Opt
% %     N = P * Ne_array(k) +1;
% %     coord = linspace(xMin, xMax, N);
% %     matA = make_K_R(P, k, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
% %     error_array(k) = calcError(P, k, A1Func,dux, matA,coord);
% %     DOF_array(k) = N;
% %     h_array(k) = L/N;
% % end
% % 
% % error_fig = plot (h_array, error_array);
% % 
% % P=3;
% % Ne_Opt=40;
% % Ne_array = linspace(1, Ne_Opt, Ne_Opt);
% % DOF_array = zeros(Ne_Opt, 1);
% % h_array = zeros(Ne_Opt, 1);
% % error_array = zeros(Ne_Opt, 1);
% % for k = 1 :1: Ne_Opt
% %     N = P * Ne_array(k) +1;
% %     coord = linspace(xMin, xMax, N);
% %     matA = make_K_R(P, k, coord, A1Func, fFunc,type_bc_l, v0, type_bc_r, vl);
% %     error_array(k) = calcError(P, k, A1Func,dux, matA,coord);
% %     DOF_array(k) = N;
% %     h_array(k) = L/N;
% % end
% % 
% % %error_fig = plot (DOF_array, error_array);
% % error_fig = plot (h_array, error_array);
% % set(gca,'FontSize',fsz);
% % xlabel('element size','FontSize', fsz);
% % ylabel('error','FontSize', fsz);
% % %title('Evoluation of numerical error with number of elements', 'FontSize', fsz);
% % leg= legend('P=1', 'P=2', 'P=3');
% % %set(leg,'FontSize',fsz)
% % 
% % 
% % % figure;
% % % error_fig = plot (h_array, error_array);
% % % set(gca,'FontSize',fsz);
% % % xlabel('element sizes','FontSize', fsz);
% % % ylabel('error','FontSize', fsz);
% % % 
% % 
