% Resolución numérica problema elíptico en 1D con condiciones dirichlet 
% homogéneas mediante elementos finitos cuadráticos 

% Parámetros del problema a*u - k*u''=f
f = @(x) 0*x+1;
a=1;
k=0.1;

% Intervalo
c=0;
d=1;

Ne = 10; % Número de elementos
h = 0.5*(d-c)/Ne;
xi = c:h:d;

run('C:\Users\Restart\Desktop\Recursos GITI\Materiales cuarto GITI\Elementos finitos\Versiones base\PEFs Elipticos\Eliptico_1d_cuad_marco\ensamblaje_mat_masrig_1d_cuad.m')
A = a*M + k*R; % Sistema "característico" de los problemas elípticos dirch. hom.

% Construcción A0 
A0 = A;
A0([1 end],:)=0;
A0(:,[1 end])=0;
A0(1,1)=1;
A0(end,end)=1;

% Cálculo b aproximando f = fh
fi = f(xi');
vect_b = M*fi ; 
vect_b([1 end]) = 0; 

uh = A0\vect_b;

plot(xi,uh)
hold on

%xi_pint = c:0.001:d;
%plot(xi_pint,sol(xi_pint), 'g')

hold off