function [ shape, dShape ] = shapeFunc

%{
***************************************************************************
INPUTS
***************************************************************************
None

***************************************************************************
VARIABLES
***************************************************************************
GQpoint     : Coordinates of 3 Gauss Quadrature points. Array [3*1].
GQweight    : Weights of 3 Gauss Quadrature points. Array [3*1].


***************************************************************************
OUTPUT
***************************************************************************
shape       : shape(i,k) is the i-th shape function evaluated at the
                     % k-th GQ point.


%}


GQpoint(1) = -sqrt(3/5);
GQpoint(2) = 0;
GQpoint(3) = sqrt(3/5);

GQweight(1) = 5/9;
GQweight(2) = 8/9;
GQweight(3) = 5/9;


shape = zeros(2,3);   

dshape = zeros(2,3);   

for k = 1:3
   ksi = GQpoint(k);
   shape(1,k)  = 0.5 * (1 - ksi);
   shape(2,k)  = 0.5 * (1 + ksi);
   dShape(1,k) = -0.5;
   dShape(2,k) = 0.5;
end

end

