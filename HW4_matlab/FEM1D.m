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
fFunc  = @(x) -256*sin(3/4*pi*x)*cos(16*pi*x);
% fFunc = @(x) 90* pi^2 * sin(3*pi*x) * sin(36*pi*x^3) -...
%     (10*sin(3*pi*x)+5)*( 216*pi*x*cos(36*pi*x^3) -...
%     11664*pi^2*x^4*sin(36*pi*x^3)) - 6480*pi^2*x^2*cos(3*pi*x)*cos(36*pi*x^3);
L=1;
xMin   =  0;
xMax   =  L;

type_bc_l=1;
v0     = 0;    % boundary value at xMin
type_bc_r=2;
vl     = 1;    % bourdary value at xMax
% type_bc_r=1;
% vl     = 0;    % bourdary value at xMax

P=1;
Tol = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adaptive meshing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ne =16;
N = P*Ne +1;
%node coordinates table
node_x = linspace(xMin, xMax, N);
%2 node numbers for each element
connectivity = zeros(Ne, 2);
for e = 1:Ne
    connectivity(e, :) = [e, e+1];
end

[K, R]= make_K_R_adp(connectivity, P, Ne, node_x, fFunc,type_bc_l, v0, type_bc_r, vl);
matA = K\R;
plot(node_x, matA, 'o');
hold on;

%---loop through to calculate element error
changed_element = true;

while changed_element == true && error > Tol
    error = 0;
    EI_tol = Tol;
    [K, R]= make_K_R_adp(connectivity, P, Ne, node_x, fFunc,type_bc_l, v0, type_bc_r, vl);
    matA = K\R;
    table_to_refine = zeros(1,1);
    for e = 1:Ne
        ce = connectivity(e,:);

        matA_elem = [matA(ce(1)); matA(ce(2))];
        coord_elem = [node_x(ce(1)), node_x(ce(2))];
        EI = elementError(matA_elem, coord_elem, xMin, xMax)
        error = error + EI*(coord_elem(2) - coord_elem(1));
        if EI > EI_tol
            % collect all the elements that need to be refined
            table_to_refine(end+1) = e;
            %---refine the mesh
            node_x(end+1) =( node_x(ce(1)) + node_x(ce(2)))/2.0;
            connectivity(e,:) = [ce(1), length(node_x)];
            connectivity(end+1, :) = [length(node_x), ce(2)];
            Ne = Ne +1;
            changed_element=true;
        end
    end
    display(table_to_refine);    
    display(error);
end
display(node_x);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot analytical solution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(node_x);
ux_arr =zeros(N, 1);
    for i=1:N
        ux_arr(i) = ana_sol_calculator(node_x(i));
    end
    
plot(node_x, ux_arr);
% legend('Ne=100', 'Ne=1000', 'Ne=10000', 'analytical solution');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %find the optimal number of element to achieve the error criteria
% % no adaptive meshing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorTol = 0.050;
% error =1;
% Ne=2;
% while error > errorTol 
%     Ne = 2* Ne;
%     N = P*Ne +1;
%     coord = linspace(xMin, xMax, N);
%     [K, R]= make_K_R(P, Ne, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%     matA = K\R;
%     error = calcError(P, Ne, matA,coord);
% end
% Ne_max = Ne;
% display(Ne_max);
% display(error);
% Ne_min = Ne/2;
% 
% while Ne_max-Ne_min > 2
%     Ne_Opt = ceil((Ne_max+Ne_min)/2);
%     N_Opt = P*Ne_Opt +1;
%     coord = linspace(xMin, xMax, N_Opt);
%     [K, R] = make_K_R(P, Ne_Opt, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%     matA = K\R;
%     error = calcError(P, Ne_Opt, matA,coord);
%     if error > errorTol
%         Ne_min = Ne_Opt;
%         %display(Ne_min);
%     else
%         Ne_max = Ne_Opt;
%         %display(Ne_max);
%     end
% end
% display(Ne_Opt);
% display(N_Opt);
% display(error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve problem 2 with a mesh pt at 1/3.0
% no adaptive meshing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N1 = 2;
% coordl = linspace(0, 1/3.0, N1);
% N2 = 4;
% coordr = linspace(1/3.0, 1, N2);
% coord = [coordl(1:length(coordl)-1), coordr];
% errorTol = 0.050;
% error =1;
% 
% while error > errorTol 
%    N1 = 2*N1;
%    coordl = linspace(0, 1/3.0, N1);
%    N2 = 2*N2;
%    coordr = linspace(1/3.0, 1, N2);
%    coord = [coordl(1:length(coordl)-1), coordr];
%    Ne = N1+N2-2;
%     [K, R]= make_K_R(P, Ne, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%     matA = K\R;
%     error = calcError(P, Ne, matA,coord);
% end
% Ne_max = Ne;
% N1_max = N1;
% N2_max = N2;
% display(Ne_max);
% display(error);
% 
% N1_min = N1/2;
% N2_min = N2/2;
% Ne_min =  N1_min+N2_min-2;
% 
% while N1_max-N1_min > 2
%    N1_Opt = ceil((N1_max+N1_min)/2);
%    coordl = linspace(0, 1/3.0, N1_Opt);
%    N2_Opt = ceil((N2_max+N2_min)/2);
%    coordr = linspace(1/3.0, 1, N2_Opt);
%    coord = [coordl(1:length(coordl)-1), coordr];
%    Ne_Opt = N1_Opt + N2_Opt -2; 
%    [K, R] = make_K_R(P, Ne_Opt, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%    matA = K\R;
%    error = calcError(P, Ne_Opt, matA,coord);
%     if error > errorTol
%         N1_min = N1_Opt;
%         N2_min = N2_Opt;
%         %display(Ne_min);
%     else
%         N1_max = N1_Opt;
%         N2_max = N2_Opt;
%         %display(Ne_max);
%     end
% end
% display(Ne_Opt);
% display(N1_Opt);
% display(N2_Opt);
% display(error);
% ************************************************************************
% Output the results and plot the solution
% ************************************************************************
% N1 = 2;
% coordl = linspace(0, 1/3.0, N1);
% N2 = 4;
% coordr = linspace(1/3.0, 1, N2);
% coord = [coordl(1:length(coordl)-1), coordr];
% 
% 
% error = zeros(6, 1); 
% for k = 2 : 10 : 50
%    N1_plot = N1 * k;
%    coordl = linspace(0, 1/3.0, N1_plot);
%    N2_plot = N2 * k;
%    coordr = linspace(1/3.0, 1, N2_plot);
%    coord = [coordl(1:length(coordl)-1), coordr];
%    Ne_plot = N1_plot + N2_plot -2; 
%    [K, R] = make_K_R(P, Ne_plot, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%    matA = K\R;
%    error(ceil(k/10)) = calcError(P, Ne_plot, matA,coord);
%     plot(coord, matA);
%     hold all on;
%     N_plot = N1_plot + N2_plot -1;
% display(N_plot);
% end
% display(error);
% 
% % Plot the numerical solution, N_ksi =10 points per element
% 
% N_plot = N1_plot + N2_plot -1;
% display(N_plot);
% u_real_array = zeros(N_plot,1);
% for k = 1: N_plot
%     u_real_array(k) = getOutput(@ana_sol_calculator,1,coord(k));
% end
% optimized_plot;
% plot(coord, u_real_array, 'r');
% fsz=30;
% xlabel('x','FontSize', fsz);
% ylabel('u','FontSize', fsz);
% leg= legend('N=11', 'N=71', 'N=131', 'N=191', 'N=251', 'analytical');
% set(leg,'FontSize',fsz)

% % %realU = ux(coord);
% % %plot (coord, realU); 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot solution for problem 2 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ne=10;
% 
%     
%     N = P*Ne +1;
%     coord = linspace(xMin, xMax, N);
%     [K, R]= make_K_R(P, Ne, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%     matA = K\R;
%     error = calcError(P, Ne, matA,coord);
% 
% 
% 
% error = zeros(6, 1); 
% for Ne = 10 : 60 : 250
%    
%     N = P*Ne +1;
%     display(N);
%    coord = linspace(xMin, xMax, N);
%     [K, R]= make_K_R(P, Ne, coord, fFunc,type_bc_l, v0, type_bc_r, vl);
%     matA = K\R;
%    error(ceil(Ne/60)) = calcError(P, Ne, matA,coord);
%     plot(coord, matA);
%     hold all on;
%     
% end
% display(error);

% Plot the numerical solution, N_ksi =10 points per element
% 
% N_plot = N;
% u_real_array = zeros(N_plot,1);
% for k = 1: N_plot
%     u_real_array(k) = getOutput(@ana_sol_calculator,1,coord(k));
% end
% optimized_plot;
% plot(coord, u_real_array, 'r');
% fsz=30;
% xlabel('x','FontSize', fsz);
% ylabel('u','FontSize', fsz);
% leg= legend('N=11', 'N=71', 'N=131', 'N=191', 'N=251', 'analytical');
% set(leg,'FontSize',fsz)

% %realU = ux(coord);
% %plot (coord, realU); 
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
