function A = make_K_sparse_fast(P, Ne, nodeCoord, A1Func, fFunc,l_BC_type, value_Left, r_BC_type, value_Right)
    
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

[GQpoint, GQweight]=lgwt(3, -1, 1);

if P==1
   N = Ne+1;
   R   = zeros(N, 1);
   II  = zeros(4*Ne,1);
   JJ  = zeros(4*Ne,1);
   VV  = zeros(4*Ne,1);
   nval =1;

   shape = {@(ksi) 0.5 * (1 - ksi);
            @(ksi) 0.5 * (1 + ksi)};
   dShape = {@(ksi) -0.5;
              @(ksi) 0.5};
end

if P==2
   N = 2*Ne+1;
    R   = zeros(N, 1);
   II  = zeros(4*Ne,1);
   JJ  = zeros(4*Ne,1);
   VV  = zeros(4*Ne,1);
   nval =1;

    %define shape functions, and ds/d(ksi)
    shape = {@(ksi) 1/2.0*(1-ksi)*ksi;
             @(ksi) (1+ksi)*(1-ksi);
             @(ksi) 1/2.0*(ksi+1)*ksi};
    dshape = {@(ksi) 1/2.0*(1-2*ksi);
              @(ksi) -2*ksi;
              @(ksi) 1/2.0*(2*ksi+1)};
end

for e=1:Ne
   x1e = nodeCoord(e);
   x2e = nodeCoord(e+1);

   % Intitialize Ke and Re to zero.
    Ke = zeros(2,2);
    Re = zeros(2,1);

    for q = 1:3   % Gauss Quadrature loop
        
      ksi = GQpoint(q);   
       % Evaluate functions of the DE at this x value.

      x = ksiTox(ksi, x1e, x2e, 1);
      a = A1Func(x); 
      f = fFunc(x);
      % Calculate Jacobian of the element
      Jacob = 0.5 * (x2e - x1e);      
      % Calculate Ke[2*2].
      for i=1:2   
          for j = 1:2
                Ke(i,j) = Ke(i,j) + (a * dShape{i}(ksi)/Jacob * dShape{i}(ksi)/Jacob) * Jacob * GQweight(q);                
          end
      end
      

       % Calculate Re        
      for i = 1:2
         Re(i) = Re(i) + f * shape{i}(ksi) * Jacob * GQweight(q);
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

%compute R

% ************************************************************************
% Apply BC at both end
% ************************************************************************
% At x = xMin
if l_BC_type==1
    R(1) = value_Left;       
    VV(1)=1 ;  %K[1,1]
    VV(3)=0 ;%K[1,2]
elseif l_BC_type ==2
     R(1) = R(1)+ value_Left;
    end
    
% At x = xMax
if r_BC_type ==1
    R(N) = value_Right; 
    VV(nval-3)= 1;
    VV(nval-1) = 0;
elseif r_BC_type ==2
    R(N) = R(N) + value_Right;
    end
        
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