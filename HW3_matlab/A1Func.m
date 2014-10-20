function A1 = A1Func(x) 
if x<0 || x>1
    error('A1Func is not defined for the given x value. ');
end
A1_array = [2.0, 2.5, 1.25, 0.25, 4.0, 1.75, 0.5, 0.75, 3.25, 1.0];
%A1_array=linspace(0.2, 0.2, 10);
A1= A1_array(max(ceil(x/0.1), 1));
end
