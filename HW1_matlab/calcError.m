function error = calcError( Ne, A1Func, AMat,k,coord)
%calculate the error of approximation  

%initialize error to 0
error =0;
a=eval(A1Func);

numeritor=0;
denominator=0;


for e=1:Ne

    du_uN_square_func = @(x) a*(-5*k/(pi)* cos(pi * k *x) + 5*x.^2 - 2/3.0 - (AMat(e+1) - AMat(e))/(coord(e+1)-coord(e))).^2;
    du_square_func = @(x) a*(-5*k/(pi)* cos(pi * k *x) + 5*x.^2 - 2/3.0).^2;
    
    u_uN_square=integral(du_uN_square_func, coord(e), coord(e+1));
    u_square = integral(du_square_func, coord(e), coord(e+1));

    numeritor = numeritor + u_uN_square;
    denominator = denominator + u_square;

end

error = sqrt(numeritor)/sqrt(denominator);

end

