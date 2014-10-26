function [dnux, unx, x_ksi] = postprocessing (ksi, AMat_element, coord_element, P)

%compute numerical solution du, u, x_ksi for any given ksi 
%on one element defined by coord_element
[shape, dshape, x_ksi] = shapeFunc(ksi, coord_element, P);

%numerical solution ux
unx = shape * AMat_element;%vectors dot product
 

%numerical solution dun/dx
Jacob = 0.5 * (coord_element(P+1) - coord_element(1));
dnux = dshape * AMat_element/Jacob;