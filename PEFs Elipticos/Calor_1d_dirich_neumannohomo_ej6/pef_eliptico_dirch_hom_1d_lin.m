% Resolución numérica problema elíptico en 1D con condiciones dirichlet 
% homogéneas y neuman no homogeneas mediante elementos finitos lineales 

% Parámetros del problema a*u - k*u''=f
f = @(x) 0*x+1;
a=1;
k=1;

% Intervalo
c=0;
d=1;

Ne = 50; % Número de elementos
h = (d-c)/Ne;
xi = c:h:d;

ensamblaje_mat_masrig_1d_lin


fi = f(xi');

B = 0*fi;
B(end)= 1; 

A = a*M + k*R; 
vect_b = M*fi+B ; 

vect_b(1) = 0; 

% Construcción A0 
A0 = A;
A0(1,:)=0;
A0(1,1)=1;

uh = A0\vect_b;

plot(xi,uh)
hold on

%xi_pint = c:0.001:d;
%plot(xi_pint,sol(xi_pint), 'g')

hold off