% Ensamblaje matrices de masas y rigideces para problema elíptico en 1D con 
% condiciones dirichlet homogéneas mediante elementos finitos lineales 

mat_masrig_elemref1d_cuad

M = sparse(2*Ne+1,2*Ne+1); %importante cuando esto crece 
R = M;

impares = 1:2:(2*Ne-1);
for i=impares
    hi = xi(i+2)-xi(i);
    M(i:i+2,i:i+2) = M(i:i+2,i:i+2) + Mg * hi*0.5;
    R(i:i+2,i:i+2) = R(i:i+2,i:i+2) + Rg * 2/hi;
end    

