% --- Script marco problemas elípticos en 2d ---

% De la formualción variacional, A(u,v), L(v), (con a, k, f, b, y las
% condiciones adicionales) definen el problema.

% Parámetros de la EDP
f = @(x,y) 2 + 0*x;
a=0;
k=1;
% b1=;
% b2=;

% Cargar coordenadas nodos
xi = p(1,:); 
yi = p(2,:);

% Matriz de conexiones filtrada y transpuesta
elem = t(1:3,:)'; 

% Detección frontera, devuelve vector de indices de nodos que forman 
% frontera, sin repeticiones.
fron_tot = unique(e(1:2,:))';

% Separar frontera en zonas según condiciones
fron_d = find(abs(0.25*xi.^2+0.0625*yi.^2 - 1) <= 0.001);
fron_n = setdiff(fron_tot,fron_d); % Diferencia de conjuntos

% Construcción fi, gi auxiliar para cond. dirichlet
fi = f(xi,yi)';
gi = 0*xi'; % Por simplicidad cero en nodos interiores
% Establecer valores en frontera según condiciones 

% Construcción matrices de rigideces, masas y convección (descomentar solo uno)
run('C:\Users\Restart\Desktop\Recursos GITI\Materiales cuarto GITI\Elementos finitos\Versiones base\PEFs Elipticos\Eliptico_2d_marco\matrices_rig_mas_conv_2d.m');
% [R M]=assema(p,t,1,1,0); % Si no se necesitan mat. convección, mucho más rápido assema

% Construcción sistema lineal característico
A = R;  
vect_b = M*fi;  
% De la formulación varaicional del prob. de elementos finitos, A y vect_b
% serán función de a, k, M, R, Cx y Cy. 

% Modificación de A y vect_b para establcer condiciones 
A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end    
vect_b(fron_d) = 0;

% Para un nodo de frontera dirichlet i, wi deja de intervenir en todas las 
% filas del sistema salvo la i, para que iguale a vect_b(i)=0. Poner a cero 
% su fila/columna salvo el termino diagonal a 1

% Resolución
wh = A0\vect_b;
trisurf(elem,xi,yi,wh)

% Comprobación frontera
hold off
figure
hold on
scatter(xi(fron_tot),yi(fron_tot),'b')
scatter(xi(fron_d),yi(fron_d),'r','.')
hold off
% Es común que la detección de los distintos trozos de frontera falle.
% Comprobarla representando. 

% --- Información extra --- 

% PDETool caracteriza el mallado mediante las matrices p (puntos), 
% e (edges), y t (conexiones). p indexa arbitrariamente los nodos y 
% lleva sus coordenadas, t y e indexan arbitrariamente triángulos y fronteras, 
% cargando los índices de los nodos que los definen


% --- Problemas frecuentes --- 

% Errores del tipo "matrix singular to working precision" pueden ser por
% vectores filas cuando corresponde columna y viceversa

% find() suele dar problemas en la comparación, poner tolerancia