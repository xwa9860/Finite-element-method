function x = shape(ksi, P, x1, x2)

%{
linear, 2nd order, 3nd order shape func
OUTPUTS:

***************************************************************************
INPUTS
***************************************************************************
p: shape function order
q: nb of guassien quadrature points used for integration

***************************************************************************
VARIABLES
***************************************************************************
GQpoint, GQweight   : Coordinates of Gauss Quadrature points. Array [Q*1]
                    : Weights of Gauss Quadrature points. Array [Q*1].


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

