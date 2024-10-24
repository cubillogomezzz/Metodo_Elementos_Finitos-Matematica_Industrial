% --- Resolución numérica MEF ec. calor 2d con condiciones dirichlet (ej.10) ---

% Parámetros del problema

f = @(x,y) 4000*exp(-10*(x.^2+y.^2));
a=0;
k=1;


xi = p(1,:); % coordenadas de los nodos
yi = p(2,:);

elem = t(1:3,:)'; % matriz de conexiones filtrada y transpuesta

fron_tot = unique(e(1:2,:))'; % indices de nodos que forman frontera, sin repeticiones

% separar frontera en zonas según condiciones 
fron_1 = find(yi==0);
fron_2 = find(xi==0);
fron_3 = setdiff(fron_tot,[fron_1 fron_2]); %diferencia de conjuntos

fi = f(xi,yi)';
% construcción de g auxiliar para cond no homogéneas:
gi = 0*xi' ; %puede valer cualquier cosa en nodos interiores, pues que valga cero
gi(fron_1)=12*xi(fron_1); %establecer valores en frontera según condiciones 
gi(fron_3)=4*xi(fron_3).^2;
%para ind_2 ya era cero

% generación matrices masas 
[R M]=assema(p,t,1,1,0);
% los 1,1,0 son a, k, 0 que matlab las pone dentro, mejor fuera, el cero ya veremos
%proceso laborioso, sencillo

A = a*M + k*R; % característico del problema
vect_b = M*fi - A*gi; 

% modificar A y vect_b estableciendo condiciones:
% para un nodo de frontera i, wi deja de intervenir en todas las filas del 
% sistema salvo la i, para que iguale a vect_b(i)=0
% poner a cero su fila/columna salvo el termino diagonal a 1

A0 = A;
A0(fron_tot,:)=0;
A0(:,fron_tot)=0;
for i=fron_tot
    A0(i,i)=1;
end    
vect_b(fron_tot) = 0;

%resolver
wh = A0\vect_b;
uh = wh + gi;

trisurf(elem,xi,yi,uh)

% comprobacion frontera
hold off
figure
hold on
scatter(xi(fron_tot),yi(fron_tot),'b')
scatter(xi(fron_1),yi(fron_1),'r','.')
scatter(xi(fron_2),yi(fron_2),'g','.')
scatter(xi(fron_3),yi(fron_3),'y','.')

hold off