function A1 = A1Func(x) 
    if x<0 || x>1
        error('A1Func is not defined for the given x value. ');
    end

    if x < 1/3.0
        A1 = 0.2;
    else
        A1 = 2.0;
    end
    
%     A1 =1.0;
    
    
end
