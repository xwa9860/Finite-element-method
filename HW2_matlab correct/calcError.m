function error = calcError( Ne, A1Func, dux, AMat, coord)
%calculate the error of approximation  
numeritor=0;
denominator=0;
for e=1:Ne
    dnux=(AMat(e+1) - AMat(e))/(coord(e+1)-coord(e));
    du_uN_square_func = @(x) A1Func(x)*(dux(x)-dnux)^2;
    du_square_func = @(x) A1Func(x)*dux(x)^2;
    
    u_uN_square=GQ_integration(1, 3, du_uN_square_func, coord(e), coord(e+1));
    u_square = GQ_integration(1, 3, du_square_func, coord(e), coord(e+1));

    numeritor = numeritor + u_uN_square;
    denominator = denominator + u_square;

end

error = sqrt(numeritor)/sqrt(denominator);

end

