function integ_value = GQ_integration(E, P, Q, func, x_min, x_max)
%{
calculate the integration value on the interval '[x_min, x_max]',
projected to some master elements [-1,1] using Gaussian Quadrature method

Input:
E: nb of elements
P: order of elements
Q: nb of GQ points
func: the function to be integrated
x_min and x_max: two ends of the integral domains

variables:
J_1D: Jacobian _ a value 1/2*[x_max-x_min] in 1D
N: nb of nodes = E+1
%}
N=E+1;
coord=linspace(x_min, x_max, N);
[GQpoint, GQweight]=lgwt(Q, -1, 1);

%shape functions for linear elements
if P==1
    %define shape functions, and ds/d(ksi)
%      shape = zeros(2,Q);   
%      dshape = zeros(2,Q);
%      for q = 1:Q
%        ksi = GQpoint(q);
%        shape(1,q)  = 0.5 * (1 - ksi);
%        shape(2,q)  = 0.5 * (1 + ksi);
%        dShape(1,q) = -0.5;
%        dShape(2,q) = 0.5;
%      end
     
     integ_value =0;
     for e=1:1:E
        x1e = coord(e);
        x2e = coord(e+1);
        J_1D = 1/2.0*(x2e-x1e);

        for q = 1:3
           ksi = GQpoint(q);
           %Inside the element e, calculate x value -- ksi of Gauss rule
           x = 0.5 * (x2e - x1e) * ksi + 0.5 * (x1e + x2e);
           integ_value = integ_value + GQweight(q)*func(x)*J_1D;
        end
     end %for e
end %end if p==1

    
    x_ksi = @(ksi) xe1*shape1(ksi) + xe2*shape2(ksi) + xe3*shape3(ksi);
end
