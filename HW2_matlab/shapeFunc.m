function [shape, dshape, x_ksi] = shapeFunc(ksi, coord_elem, P)
%{
INPUT:
coord_elem is the P+1 coordinates in one elements
P:order of basis functions 

OUTPUT:
shape: array of values [phi_hat1(ksi), phi_hat2(ksi)...] for a given ksi
dshape: array of values [dphi_hat1(ksi), dphi_hat2(ksi)...] for a given ksi
x_ksi: the value of x for an corresponding ksi in the element defined by
coord_elem
%}

if nargin ==3
else
    error('need 3 input arguments, got %d', nargin)
end
   
[m,n] = size(coord_elem);
if n== P+1
else
    error('coord_elem size %d given to shapeFunc does not match order %d, should give P+1 coordinates for an element', n, P);
end

%define shape functions(phi_hat), and their derivatives ds/d(ksi)
if P==1
   shapeFunc = {@(ksi) 0.5 * (1 - ksi);
            @(ksi) 0.5 * (1 + ksi)};
   dshapeFunc = {@(ksi) -0.5;
              @(ksi) 0.5};
          
elseif P==2
    
    shapeFunc = {@(ksi) -1/2.0*(1-ksi)*ksi;
             @(ksi) (1+ksi)*(1-ksi);
             @(ksi) 1/2.0*(ksi+1)*ksi};
    dshapeFunc = {@(ksi) -1/2.0*(1-2*ksi);
              @(ksi) -2*ksi;
              @(ksi) 1/2.0*(2*ksi+1)};
elseif P==3
    shapeFunc = { @(ksi) -1/16.0*(ksi-1)*(3*ksi+1)*(3*ksi-1);
            @(ksi) 9/16.0*(ksi-1)*(ksi+1)*(3*ksi-1);
            @(ksi) -9/16.0*(ksi-1)*(ksi+1)*(ksi*3+1);
            @(ksi) 1/16.0*(ksi+1)*(3*ksi+1)*(3*ksi-1)
        } ;
    dshapeFunc = {@(ksi) -1/16.0*(27*ksi^2-18*ksi-1);
        @(ksi) 9/16.0*(9*ksi^2-2*ksi-3);
        @(ksi) -9/16.0*(9*ksi^2+2*ksi-3);
        @(ksi) 1/16.0*(27*ksi^2+18*ksi-1)
        };
else
    error('this function only support P =1,2 or 3');
end

x_ksi=0;

shape = zeros(1, P+1);
dshape = zeros(1, P+1);
for p=1:1:P+1
    shape(p) = shapeFunc{p}(ksi); 
    dshape(p) = dshapeFunc{p}(ksi);
    x_ksi = x_ksi+ coord_elem(p)*shapeFunc{p}(ksi);
end


end