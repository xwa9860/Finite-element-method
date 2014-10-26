function [ux, dux] = ana_sol_calculator(x)
%{
define the analytical solution
%}

if x<0 || x>1
    error('ux and dux is not defined for the given x value. only defined on x=[0,1]');
end

seg_nb =2;
C1_array = zeros(seg_nb, 1); 
C1_array(1) = 5 + 7680*sqrt(2) /4087 / pi;
C1_array(2) = 0.5 + 768*sqrt(2) /4087 / pi;

C2_array = zeros(seg_nb, 1);
C2_array(1) =1/0.2*512* (268*sin(61*pi*x/4)/61/pi -244*sin(67*pi*x/4)/67/pi)/4087/pi;
C2_array(2) =1/2.0*512* (268*sin(61*pi*x/4)/61/pi -244*sin(67*pi*x/4)/67/pi)/4087/pi+...
2304*sqrt(2)*(8210 - 768*sqrt(3)+4087*pi)/16703569 /pi^2 + 1.5;

if x < 1/3.0 
    ux = C1_array(1)*x+ C2_array(1);
    dux = C1_array(1) + 1/0.2*512* (268/4*cos(61*pi*x/4) -244/4*cos(67*pi*x/4))/4087/pi;
else
    ux = C1_array(2)*x+ C2_array(2);
    dux = C1_array(2) + 1/2.0*512* (268/4*cos(61*pi*x/4) -244/4*cos(67*pi*x/4))/4087/pi;
end


% ux = (10*sin(3*pi*x) + 5)*sin(36*pi*x^3);
% dux = 6*pi*(5*sin(36*pi*x^3)*cos(3*pi*x) + 90* x^2 * (2*sin(3*pi*x) +1)*cos(36*pi*x^3));

end



