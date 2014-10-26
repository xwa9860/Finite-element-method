function EI = elementError (AMat_elem, coord_elem, xMin, xMax)
P=1;
Q=5;
coord_elem_ksi = [-1.0, 1.0]; 
hI = coord_elem(2) - coord_elem(1);
L = xMax-xMin;    

    A1= @(ksi)A1Func(getOutput(@shapeFunc, 3, ksi, coord_elem, P));
    dux_true =@(ksi) getOutput(@ana_sol_calculator, 2, (getOutput(@shapeFunc, 3, ksi, coord_elem, P)));
    dux_N = @(ksi) getOutput(@postprocessing, 1, ksi, AMat_elem, coord_elem, P);
    
    J = hI/2;
    du_uN_square_func = @(ksi) A1(ksi)*(dux_true(ksi)-dux_N(ksi))^2*J;
    u_uN_square=GQ_integration(P, Q, du_uN_square_func, coord_elem_ksi);
    numeritor =  u_uN_square;
    
    %--- calculate demoninator
    denominator=0;
    Ne =1000;                       % only for calculate the demoninator
    N = Ne+1;
    coord = linspace(xMin, xMax, N);% doesn't depend on the mesh used for solving the equation
    for e=1:Ne
        Nl= P*(e-1)+1;
        Nr= P*e+1;
        coord_elem= coord(Nl:Nr);
         J = (coord(Nr)-coord(Nl))/2;

        x = @(ksi) getOutput(@shapeFunc, 3, ksi, coord_elem, P);
        A1= @(ksi)A1Func(x(ksi));
        dux_true =@(ksi) getOutput(@ana_sol_calculator, 2, x(ksi));
        du_square_func = @(ksi) A1(ksi)*(dux_true(ksi))^2*J;
        u_square = GQ_integration(P, Q, du_square_func, coord_elem_ksi);

        denominator = denominator + u_square;

    end

    EI = (numeritor/hI)/(denominator/L);



end
