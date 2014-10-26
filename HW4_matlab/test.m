%analytical solution for the problem
% x=1;
% % C1_array = [2.0, 2.5, 1.25, 0.25, 4.0, 1.75, 0.5, 0.75, 3.25, 1.0];
% % C1= C1_array(max(ceil(x/0.1), 1));
% % 
% % C2_array = [];
% % C2= C2_array(max(ceil(x/0.1), 1));
% 
x_array= linspace(0,1,11);

dux_func =@(x, C1) 1/A1Func(x)*(512*(67*cos(61*pi*x/4) -61*cos(67*pi*x/4))/4087/pi)+C1/A1Func(x);
C1_array = zeros(10, 1);

%initial value for du/dx at n=11 : x=1 is A1dux =1
x_11 =1;
A1dux=1*A1Func(x_11);

for n=11:-1:2

    x_n = 0.1*(n-1); %coordinates at Nth point N=1:10
    syms C1
    s= solve( A1Func(x_n)*dux_func(x_n, C1) == A1dux, C1);
    C1_array(n-1) = s;
    
    %go to the left end of the segment
    x_n_1 = 0.1*(n-2); % x(n-1)
    A1dux = dux_func(x_n_1, s)*A1Func(x_n);%interface x(n-1)calculate dux with C1 value at this segment
end
display(C1_array);

C2_array = zeros(10, 1);
ux_func =@(x, C1, C2) -1/A1Func(x)*(137216*sin(61*pi*x/4)/61/pi -124928*sin(67*pi*x/4)/67/pi)/4087/pi+...
      C1*x/A1Func(x) + C2;

ux=0;

for n =1:1:10
  syms C2
  s = solve(ux_func(x_array(n), C1_array(n), C2) == ux, C2);
  C2_array(n) = s;
    
  %go to the right end of the segment
  ux = ux_func(x_array(n+1), C1_array(n), s)*A1Func(x_array(n));%interface x(n-1)calculate dux with C1 value at this segment
end
display(C2_array);


% x = linspace(0, 1, 100);
% plot(x, u1(x));

% hold on;
% ux =@(x) 5*((4087 * pi *x + 1536*sqrt(2) *x + ...
%     137216*sin(61*pi*x/4)/61/pi -124928*sin(67*pi*x/4)/67/pi)...
%     /4087/pi);


plot(x_arr, ux(x));
% syms u(t)
% Du(t) = diff(u);
% u1(t) = dsolve(diff(u,t,t) == 256*sin(3/4*pi*t)*cos(16*pi*t), u(0) == 0, Du(0) == 0, t)


