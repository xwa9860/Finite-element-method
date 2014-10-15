function integ_value = GQ_integration(P, Q, func, coord_element)
%{
calculate the integration value of  function 'func'on the interval '[coord_element]',
by projecting it to some master elements [-1,1] using Gaussian Quadrature method

LIMIT:
the function should be approximated by the sum of the integral on E elements

Input:
P: order of elements
Q: nb of GQ points
func: the function to be integrated
x_min and x_max: two ends of the integral domains

variables:
J_1D: Jacobian _ a value 1/2*[x_max-x_min] in 1D

%}
[GQpoint, GQweight]=lgwt(Q, -1, 1);
     
integ_value =0;
J_1D = 1/2.0*(coord_element(1+P)-coord_element(1));

for q = 1:Q
    ksi = GQpoint(q);
    %Inside the element e, calculate x value -- ksi of Gauss rule
    [shape, dshape, x] = shapeFunc(ksi, coord_element, P);
    integ_value = integ_value + GQweight(q)*func(x)*J_1D;
end
    
end 

