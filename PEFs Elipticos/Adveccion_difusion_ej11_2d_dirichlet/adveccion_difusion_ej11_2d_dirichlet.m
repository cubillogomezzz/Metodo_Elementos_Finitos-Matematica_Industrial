% --- Resolución numérica MEF ec. adveccion-difusión (tipo ej.11) 2d con condiciones dirichlet ---

load('entorno_mallado.mat')


f = @(x,y) 0*x;
a=0;
k=1;
b1 = 2;
b2 = 1;

xi = p(1,:); 
yi = p(2,:);
elem = t(1:3,:)';
fron_tot = unique([e(1,:) e(2,:)]); 
tol = 0.001;
fron_1 = find((yi -(0.5 - xi))<=tol);



fi = f(xi,yi)';
gi = 0*xi';
gi(fron_1)=1;

% sacar matrices de rigideces, masas y convección
run('C:\Users\Restart\Desktop\Recursos GITI\Materiales cuarto GITI\Elementos finitos\Versiones base\PEFs Elipticos\Eliptico_2d_marco\matrices_rig_mas_conv_2d.m');

A = b1*Cx + b2*Cy+k*R;
vect_b = -A*gi; 

A0 = A;
A0(fron_tot,:)=0;
A0(:,fron_tot)=0;
for i=fron_tot
    A0(i,i)=1;
end    
vect_b(fron_tot) = 0;

wh = A0\vect_b;
uh = wh + gi; 


trisurf(elem,xi,yi,uh)
%shading interp

% comprobacion frontera
hold off
figure
hold on
scatter(xi(fron_tot),yi(fron_tot),'b')
scatter(xi(fron_1),yi(fron_1),'r','.')
hold off