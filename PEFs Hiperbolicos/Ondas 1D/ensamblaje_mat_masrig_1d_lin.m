% Ensamblaje matrices de masas y rigideces para problema elíptico en 1D con 
% condiciones dirichlet homogéneas mediante elementos finitos lineales 

mat_masrig_elemref1d_lin

M = sparse(Ne+1,Ne+1); %importante cuando esto crece 
R = M;

for i=1:Ne

    hi = xi(i+1)-xi(i);
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + Mg * hi*0.5;
    R(i:i+1,i:i+1) = R(i:i+1,i:i+1) + Rg * 2/hi;

end    

