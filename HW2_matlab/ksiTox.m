function x = ksiTox(ksi, x1, x2, P)

if P==1
 %Inside the element e, calculate x value -- ksi of Gauss rule
      x = 0.5 * (x2 - x1) * ksi + 0.5 * (x1 + x2);
end

end