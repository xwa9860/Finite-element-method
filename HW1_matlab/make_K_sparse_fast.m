function A = make_K_sparse_fast(N, k,nodeCoord, A1Func, fFunc,uLeft,uRight)
    
%{
***************************************************************************
INPUTS
***************************************************************************
N           :number of nodes

***************************************************************************
VARIABLES
***************************************************************************
GQpoint     : Coordinates of 3 Gauss Quadrature points. Array [3*1].
GQweight    : Weights of 3 Gauss Quadrature points. Array [3*1].
Ke          : Elemental stiffness matrix [2x2].
Fe          : Elemental load vector [2x1].

***************************************************************************
OUTPUT
***************************************************************************
K       : stiffness matrix. Sparse matrix [N*N]
R       : load vector. Array [1*N]


%}

% ************************************************************************
% Table of  GQ points and weights
% ************************************************************************

GQpoint(1) = -sqrt(3/5);
GQpoint(2) = 0;
GQpoint(3) = sqrt(3/5);

GQweight(1) = 5/9;
GQweight(2) = 8/9;
GQweight(3) = 5/9;

% ************************************************************************
% Compute sparse matrix K and R
% ************************************************************************

% Initialization
Ne  = N-1;
R   = zeros(N, 1);
II  = zeros(4*Ne,1);
JJ  = zeros(4*Ne,1);
VV  = zeros(4*Ne,1);
nval =1;


[S, dS ] = shapeFunc;



for e=1:Ne
   x1e = nodeCoord(e);
   x2e = nodeCoord(e+1);

   % Intitialize Ke and Re to zero.
    Ke = zeros(2,2);
    Re = zeros(2,1);

    for q = 1:3   % Gauss Quadrature loop

      ksi = GQpoint(q);   
      %Inside the element e, calculate x value -- ksi of Gauss rule
      x = 0.5 * (x2e - x1e) * ksi + 0.5 * (x1e + x2e);
      % Evaluate functions of the DE at this x value.
      a = eval(A1Func);   
      f = eval(fFunc);
      % Calculate Jacobian of the element
      Jacob = 0.5 * (x2e - x1e);

      % Calculate Ke[2*2].
       for i = 1:2
         for j = 1:2
            Ke(i,j) = Ke(i,j) + (a * dS(i,q)/Jacob * dS(j,q)/Jacob) * Jacob * GQweight(q);
         end
       end

       % Calculate Re        
      for i = 1:2
         Re(i) = Re(i) + f * S(i,q) * Jacob * GQweight(q);
      end
          
    end% End of GQ loop

    R(e) = R(e)+Re(1);
    R(e+1) = R(e+1)+Re(2);
    II(nval) = e;
    JJ(nval) = e;
    VV(nval) = Ke(1,1);
    II(nval+1) = e+1;
    JJ(nval+1) = e+1;
    VV(nval+1) = Ke(2,2);
    II(nval+2) = e;
    JJ(nval+2) = e+1;
    VV(nval+2) = Ke(1,2);
    II(nval+3) = e+1;
    JJ(nval+3) = e;
    VV(nval+3) = Ke(2,1);
    nval = nval + 4;
    
end%all element

% ************************************************************************
% Apply diriclet BC at both end
% ************************************************************************
% At x = xMin
R(1) = uLeft;       
VV(1)=1 ;  
VV(3)=0 ;%K[1,2]

% At x = xMax
R(N) = uRight; 
display(II(nval-3));
display(JJ(nval-3));
VV(nval-3)= 1;
VV(nval-1) = 0;

% Delete unused indices:
% II = II(1:nval);
% JJ = JJ(1:nval);
% VV = VV(1:nval);

%Assemble the sparse matrix all at once.
% K(i,j)=V,If you repeat indices, e.g. say I,J,V = [ 1 1 ], [ 1 1 ], [2.0
% 3.0], the assembly will add repeated values, so you'll end up
% with a 5.0 at 1,1 in the above example.
K = sparse( II, JJ, VV, N,N);

A = K\R;     
end