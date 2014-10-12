
clc;
clear all;
close all;
%test GQ_integration function
func_f=@(x) x^2;
integ = GQ_integration(1, 3, func_f, 0, 2);
error = abs(integ - 8/3.0);
if error<0.001
    display('integration function is good');
else
    display('integration function is bad, error is ' );
    display(error);
end

%test shape functions
%ksi=[-1, -1/3.0, 1/3.0, 1];
ksi=[-1,0,1];
P=2;
p=P+1;
coordE=[1 2 3];
shape_1=zeros(1,p);
shape_2=zeros(1,p);
shape_3=zeros(1,p);
shape_4=zeros(1,p);
x_coord=zeros(1,p);
for p =1:1:P+1
[s, ds, x] = shapeFunc(ksi(p), coordE, P );
shape_1(p) = s(1);
shape_2(p) = s(2);
shape_3(p) = s(3);
%shape_4(p) = s(4);
x_coord(p)=x;
end
hold on;
plot(ksi, shape_1, 'oy');
plot(ksi, shape_2, 'or');
plot(ksi, shape_3, 'ob');
%plot(ksi, shape_4, 'og');
display(x_coord);

nodecoord = linspace(0, 1, 10);
P=1;
e=2;
%index_l=0;
index_l =e+1-1;
%index_r=0;
index_r=2*e-1+P;
coordElem = nodecoord(e:4)
coordElement = nodecoord(index_l:2*e-1+P)