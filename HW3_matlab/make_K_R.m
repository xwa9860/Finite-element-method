function [KE_table, R] = make_K_R(P, Ne, nodeCoord, fFunc, l_BC_type, value_Left, r_BC_type, value_Right)
    
%{
***************************************************************************
INPUTS
***************************************************************************
P   : polynomial order for basis function on each element,
Ne  : number of mesh elements
nodeCoord: coordinates array for all the nodes
%}
Q=5;%gaussian quadrature point number
[GQpoint, GQweight]=lgwt(Q, -1, 1);

N = P*Ne+1; %number of nodes

% ************************************************************************
% Compute sparse matrix K and R
% ************************************************************************
% Initialization
R   = zeros(N, 1);
KE_table = zeros(N, 3);
% II  = zeros((P+1)*(P+1)*Ne,1);
% JJ  = zeros((P+1)*(P+1)*Ne,1);
% VV  = zeros((P+1)*(P+1)*Ne,1);
% nval =1;
R_index=1;

for e=1:Ne
    Nl= P*(e-1)+1;
    Nr= P*e+1;
    coordElement = nodeCoord(Nl:Nr);

    % Intitialize Ke and Re to zero.
    Ke = zeros(P+1, P+1);
    Re = zeros(P+1,1);

    for q = 1:Q   % Gauss Quadrature loop

        ksi = GQpoint(q);   
        [shape, dShape, x_ksi] = shapeFunc(ksi, coordElement, P);

        % Calculate Jacobian of the element
        Jacob = 0.5 * (nodeCoord(Nr) - nodeCoord(Nl));

    %     % Calculate Ke
    %     for i=1:(P+1)   
    %         for j = 1:(P+1)
    %             Ke(i,j) = Ke(i,j) + (-A1Func(x_ksi) *dShape(i)/Jacob * dShape(j)/Jacob) * Jacob * GQweight(q);  
    %        end
    %     end
        i=1; j=1;
        KE_table (e, 1) = KE_table(e, 1)+(-A1Func(x_ksi) *dShape(i)/Jacob * dShape(j)/Jacob) * Jacob * GQweight(q);  
        i=1; j=2;
        KE_table (e, 2) = KE_table(e, 1)+(-A1Func(x_ksi) *dShape(i)/Jacob * dShape(j)/Jacob) * Jacob * GQweight(q);  
        i=2; j=2;
        KE_table (e, 3) = KE_table(e, 1)+(-A1Func(x_ksi) *dShape(i)/Jacob * dShape(j)/Jacob) * Jacob * GQweight(q);    
    end
    
    % Calculate Re
    for i = 1:(P+1)
        Re(i) = Re(i) + fFunc(x_ksi) * shape(i) * Jacob * GQweight(q);
    end

    end% End of GQ loop
  
%     for i=1:(P+1)
%         for j=1:(P+1)
%             II(nval) = (e-1)*P+i;
%             JJ(nval) = (e-1)*P+j;
%             VV(nval) = Ke(i,j);
%             nval =nval+1;
%         end
%     end
    
    for k=0:1:P
    R(R_index) = R(R_index)+Re(k+1);
    R_index=R_index+1;
    end
    R_index=R_index-1;


end%all element


% ************************************************************************
% Apply BC at both end
% ************************************************************************

% At x = xMin
if l_BC_type==1
    R(1) = value_Left;
    
    KE_table(1,2) = 0;
    KE_table(1,1) = 1;
    
%     VV(1)=1;
%     if P==1        
%         VV(2)=0 ;  %K[1,:] 
%     elseif P==2
%         VV(2)=0 ;  %K[1,:]
%         VV(3)=0;
%     elseif P==3
%         VV(2)=0;
%         VV(3)=0;
%         VV(4)=0;
%     else
%         error('does not support higher order than P=3');
%     end%if
elseif l_BC_type ==2
     R(1) = R(1) - value_Left;    
end
    
% At x = xMax
if r_BC_type ==1
    R(N) = value_Right;
    
    KE_table(N, 3)=1;
    KE_table(N, 2)=0;
%     if P==1      
%         VV(nval-2)= 0;%K[N,:]   
%     elseif P==2
%         VV(nval-2)= 0;
%         VV(nval-3)=0;
%     elseif P==3
%         VV(nval-2)= 0;
%         VV(nval-3)=0;
%         VV(nval-4)=0;
%     else
%         error('does not support P higher than 3')
%     end
%     
%     VV(nval-1) = 1; %K[N,N]
elseif r_BC_type ==2
    R(N) = R(N) - value_Right;    
end
        

% Delete unused indices:
% II = II(1:nval);
% JJ = JJ(1:nval);
% VV = VV(1:nval);

%Assemble the sparse matrix all at once.
% K(i,j)=V,If you repeat indices, e.g. say I,J,V = [ 1 1 ], [ 1 1 ], [2.0
% 3.0], the assembly will add repeated values, so you'll end up
% with a 5.0 at 1,1 in the above example.
%K = sparse( II, JJ, VV, N,N);



end