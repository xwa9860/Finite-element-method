function [ux, dux] = ana_sol_calculator(x)
%{
define the analytical solution
%}
A1_array = transpose([2.0, 2.5, 1.25, 0.25, 4.0, 1.75, 0.5, 0.75, 3.25, 1.0]);
%A1_array=transpose(linspace(0.2, 0.2, 10));
if x<0 || x>1
    error('ux and dux is not defined for the given x value. only defined on x=[0,1]');
end

x_array= linspace(0,1,11);

A1dux_func =@(x, C1) (512*(67*cos(61*pi*x/4) -61*cos(67*pi*x/4))/4087/pi)+C1;
C1_array = zeros(10, 1);

%boundary value for du/dx at n=11 : x=1 is A1dux =1
A1dux=1;

for n=11:-1:2
    syms C1
    s= solve(A1dux_func(x_array(n), C1) == A1dux, C1);
    C1_array(n-1) = s;
    
    %go to the left end of the segment
    A1dux = A1dux_func(x_array(n-1), s);%interface x(n-1)+ calculate A1dux 
end

%display(C1_array);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calc C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C2_array = zeros(10, 1);
A1ux_func =@(x, C1, C2) (137216*sin(61*pi*x/4)/61/pi -124928*sin(67*pi*x/4)/67/pi)/4087/pi+...
      C1*x+ C2;
%ux at x=0 is 0
ux=0;
for n =1:1:10
  syms C2
  s = solve(A1ux_func(x_array(n), C1_array(n), C2)/A1_array(n) == ux, C2);
  C2_array(n) = s;    
  %go to the right end of the segment
  ux = A1ux_func(x_array(n+1), C1_array(n), s)/A1_array(n);%interface x(n+1)- calculate dux with C1 value at this segment
end

%display(C2_array./A1_array);

ux =A1ux_func(x, C1_array(max(ceil(x/0.1), 1)), C2_array(max(ceil(x/0.1), 1)))/A1Func(x);

dux = A1dux_func(x, C1_array(max(ceil(x/0.1), 1)))/A1Func(x);


