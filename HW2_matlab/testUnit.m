
%test GQ_integration function
func_f=@(x) x^2;
integ = GQ_integration(100, 1, 3, func_f, 0, 2);
error = abs(integ - 8/3.0);
if error<0.001
    display('integration function is good');
else
    display('integration function is bad, error is ' );
    display(error);
end

ksi=[-1:0.2:1]
x = ksiTox(ksi, 1, 3, 1);
plot(ksi, x, 'o');
