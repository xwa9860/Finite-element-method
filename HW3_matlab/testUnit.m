
clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test GQ_integration function
func_f=@(x) x^2;
coordE= [0, 2/3.0, 4/3.0, 2];
integ = GQ_integration(3, 3, func_f, coordE);
error = abs(integ - 8/3.0);
if error<0.00001
    display('integration function is good');
else
    display('integration function is bad, error is ' );
    display(error);
end

coordE= [-1, -1/3.0, 1/3.0, 1];
integ = GQ_integration(3, 3, func_f, coordE);
error = abs(integ - 2/3.0);
if error<0.00001
    display('integration function is good');
else
    display('integration function is bad, error is ' );
    display(error);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test getOutput function
x = getOutput(@testFunc, 1, 1,2);
y = getOutput(@testFunc, 2, 1,2);
if x~=1|y~=2
    display('getOutput function is wrong')
else
    display('getOutput function is good')
end

A1Func=@(x) x;
ksi =-1;
a= A1Func(getOutput(@shapeFunc, 3, ksi, [1,2], 1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
hold off;
display(x_coord);

nodecoord = linspace(0, 1, 10);
P=1;
e=2;
%index_l=0;
index_l =e+1-1;
%index_r=0;
index_r=2*e-1+P;
coordElem = nodecoord(e:4);
coordElement = nodecoord(index_l:2*e-1+P);

%plot shape functions
N_ksi =100;
P=2;
index =1;
Phi_ksi_array = zeros(N_ksi,(P+1));
ksi_array = linspace(-1, 1, N_ksi);
coord_elem = linspace(1,2, P+1);


    
    for k = 1:1:N_ksi
        [shape, dshape, x_ksi] = shapeFunc(ksi_array(k), coord_elem, P);
        x_ksi_array(index) = x_ksi;
        Phi_ksi_array(index,:) = shape;
        index = index +1;
    end

    % Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 30;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
    
    
    for p = 1:(P+1)
        plot(x_ksi_array, Phi_ksi_array(:,p));
        hold all;
    end
    title('Quadratic shape functions')
    xlabel('x')
    ylabel('\phi(x)')
    legend('\phi_1','\phi_2','\phi_3','\phi_4');
    hold off;
set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
set(gca,'YTick',0:10); %<- Still need to manually specific tick marks

