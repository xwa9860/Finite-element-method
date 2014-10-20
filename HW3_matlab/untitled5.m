ux =@(x) 5*((4087 * pi *x + 1536*sqrt(2) *x + ...
    137216*sin(61*pi*x/4)/61/pi -124928*sin(67*pi*x/4)/67/pi)...
    /4087/pi);

x = linspace(0, 1, 1000);
u=zeros(1000,1);
for i =1:1:1000
u(i)=  getOutput(@ana_sol, 1, x(i));
end

plot(x, u);
hold all on;
plot(x, ux);

legend('u','ux');