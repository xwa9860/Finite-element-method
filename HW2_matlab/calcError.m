function error = calcError( Ne, A1Func, AMat, coord)
%calculate the error of numerical approximation  

%initialize error to 0
error =0;
a= A1Func(x);
numeritor=0;
denominator=0;


for e=1:Ne

    du_uN_square_func = @(x) a*(5*((4087 * pi  - 32768*sqrt(2)  + 61*pi/4*137216*cos(61*pi*x/4)/61/pi -124928*67*pi/4*cos(67*pi*x/4)/67/pi)/4087/pi) - (AMat(e+1) - AMat(e))/(coord(e+1)-coord(e))).^2;
    du_square_func = @(x) a*(5*((4087 * pi  - 32768*sqrt(2)  + 61*pi/4*137216*cos(61*pi*x/4)/61/pi -124928*67*pi/4*cos(67*pi*x/4)/67/pi)/4087/pi)).^2;
    
    u_uN_square=integral(du_uN_square_func, coord(e), coord(e+1));
    u_square = integral(du_square_func, coord(e), coord(e+1));

    numeritor = numeritor + u_uN_square;
    denominator = denominator + u_square;

end

error = sqrt(numeritor)/sqrt(denominator);

end

