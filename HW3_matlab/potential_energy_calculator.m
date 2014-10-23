function J = potential_energy_calculator (P, f, Ne, AMat, coord)

Q=5;

if P==1
    coord_elem_ksi = [-1.0, 1.0];
elseif P==2
    coord_elem_ksi = [-1,0,1];
elseif P==3
    coord_elem_ksi = [-1, -1/3.0, 1/3.0, 1];
end

u_square =0;
fu=0;

for e=1:Ne
    Nl= P*(e-1)+1;
    Nr= P*e+1;
    coord_elem= coord(Nl:Nr);
    AMat_elem = AMat(Nl:Nr);
   
    du_square_func = @(ksi) A1Func(getOutput(@shapeFunc, 3, ksi, coord_elem, P))*(getOutput(@ana_sol_calculator, 2, (getOutput(@shapeFunc, 3, ksi, coord_elem, P))))^2;
    fu_func = @(ksi) f(getOutput(@shapeFunc, 3, ksi, coord_elem, P))* getOutput(@ana_sol_calculator, 1, (getOutput(@shapeFunc, 3, ksi, coord_elem, P)));
    u_square = u_square + GQ_integration(P, Q, du_square_func, coord_elem_ksi);
    fu = fu + GQ_integration(P, Q, fu_func, coord_elem_ksi);
end

J = 1/2.0 * u_sqaure + fu + 1*getOutput(@ana_sol_calculator, 1, 1);

end