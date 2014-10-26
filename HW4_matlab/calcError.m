function error = calcError( P, Ne, AMat, coord)
%calculate the error of approximation  
numeritor=0;
denominator=0;

Q=5;
if P==1
    coord_elem_ksi = [-1.0, 1.0];
elseif P==2
    coord_elem_ksi = [-1,0,1];
elseif P==3
    coord_elem_ksi = [-1, -1/3.0, 1/3.0, 1];
end

for e=1:Ne
    Nl= P*(e-1)+1;
    Nr= P*e+1;
    coord_elem= coord(Nl:Nr);
    AMat_elem = AMat(Nl:Nr);
    J = (coord(Nr)-coord(Nl))/2;
    
    x = @(ksi) getOutput(@shapeFunc, 3, ksi, coord_elem, P);
    A1= @(ksi)A1Func(x(ksi));
    dux_true =@(ksi) getOutput(@ana_sol_calculator, 2, x(ksi));
    dux_N = @(ksi) getOutput(@postprocessing, 1, ksi, AMat_elem, coord_elem, P);
    
    du_uN_square_func = @(ksi) A1(ksi)*(dux_true(ksi)-dux_N(ksi))^2*J;
    du_square_func = @(ksi) A1(ksi)*(dux_true(ksi))^2*J;
    
    u_uN_square=GQ_integration(P, Q, du_uN_square_func, coord_elem_ksi);
    u_square = GQ_integration(P, Q, du_square_func, coord_elem_ksi);

    numeritor = numeritor + u_uN_square;
    denominator = denominator + u_square;

end

error = sqrt(numeritor)/sqrt(denominator);

end

