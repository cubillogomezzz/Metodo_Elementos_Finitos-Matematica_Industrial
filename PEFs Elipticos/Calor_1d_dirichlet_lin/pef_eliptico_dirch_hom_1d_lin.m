% Resolución numérica problema elíptico en 1D con condiciones dirichlet 
% homogéneas mediante elementos finitos lineales 

% Parámetros del problema a*u - k*u''=f
f = @(x) 0*x+1;
a=1;
k=0.1;

% Intervalo
c=0;
d=1;

Ne = 5; % Número de elementos
h = (d-c)/Ne;
xi = c:h:d;

ensamblaje_mat_masrig_1d_lin
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