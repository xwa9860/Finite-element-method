function A1 = A1Func(x) 
A1_array = [2.0, 2.5, 1.25, 0.25, 4.0, 1.75, 0.5, 0.75, 3.25, 1.0];
A1= A1_array(max(ceil(x/0.1), 1));
end
