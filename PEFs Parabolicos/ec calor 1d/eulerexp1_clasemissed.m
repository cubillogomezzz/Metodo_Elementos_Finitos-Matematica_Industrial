% problema parabolico ej1 hoja 6 (ec calor)
% método euler implícito (sí es imp) 

% todo cambia mucho para cada problema, ver formulacion algebraica



% Resolución numérica problema elíptico en 1D con condiciones dirichlet 
% homogéneas mediante elementos finitos lineales 

% Parámetros del problema a*u - k*u''=f
f = @(x,t) 0*x + 0*t;
a=0;
k=0.1;

dt = 0.1;
t0 = 0;
tf = 100;
Nt = (tf-t0)/dt;

% Intervalo
c=0;
d=1;

Ne = 10; % Número de elementos
h = (d-c)/Ne;
xi = c:h:d;

ensamblaje_mat_masrig_1d_lin

A = (1+a*dt)*M + k*dt*R; % Sistema "característico" de los problemas elípticos dirch. hom. y metodo euler imp

% Construcción A0 
A0 = A;
A0([1 end],:)=0;
A0(:,[1 end])=0;
A0(1,1)=1;
A0(end,end)=1;

% Cálculo b aproximando f = fh

uhn = sin(pi*x').^100; %cond inicial
for n = 1:Nt
    fi = f(xi',n*dt);
    vect_b = M*uhn + dt*M*fi; % caract del metodo y problema %probs dimension
    vect_b([1 end]) = 0; 

    uhn = A0\vect_b; %si hiciese falta guardamos cada un0

    plot(xi,uhn)
    pause
    hold on
end



%xi_pint = c:0.001:d;
%plot(xi_pint,sol(xi_pint), 'g')

hold off