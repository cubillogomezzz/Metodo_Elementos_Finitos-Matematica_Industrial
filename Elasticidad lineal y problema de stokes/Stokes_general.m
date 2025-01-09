% --- Problema de Stokes con dirichlet general, neuman flujo libre --- 

% Para satisfacer condicion Inf-Sup de existencia y unicidad, Vh cuadrático
% para velocidad, lineal para presión. 

nu = 1;
rho = 1;

f1 = @(x,y) 0*x; 
f2 = @(x,y) 0*x;

if 0 %para no estar ejecutando varias veces
    
    myPDE = createpde;  
    
    % Generamos geometria en pdetool, hacemos boundary -> export decomposed geometry 
    geometryFromEdges(myPDE,g); 
    % Parametros mallado cuadrático, hmax longitud característica máxima
    generateMesh(myPDE,'hmax',0.1,'geometricorder','quadratic'); 
    myPDE.Mesh
    
    xi = myPDE.Mesh.Nodes(1,:);
    yi = myPDE.Mesh.Nodes(2,:);
    
    elem = myPDE.Mesh.Elements(1:3,:)'; % Cuidado transponer
    
    applyBoundaryCondition(myPDE,'neumann','Edge',1:size(g,2),'q',1) 
    % Establecemos todo Neuman aunque no lo sea para el truco de detectar
    % frontera. Impondremos condiciones manualmente.
    specifyCoefficients(myPDE,'a',1,'c',1,'m',0,'d',0,'f',0); 
    % Coeficientes neutros para tratar desde fuera. En nuestros ejercicios
    % son todos 1. 


    FEM_KM = assembleFEMatrices(myPDE,'none'); % 'none' evita ensablar condiciones en matrices
    
    M = FEM_KM.A;
    R = FEM_KM.K;
    
    [Cxl ,Cyl] = calcular_matrices_stokes(myPDE); % Solo soporta geometrías poligonales
    
end


fron_tot = find(diag(FEM_KM.Q)>0);
% Matriz de masas sobre la frontera, si vemos su diagonal donde sea positiva, 
% nos devuelve nodos de la frontera. Cuidado el size que saca esto, transponer si es necesario 

fron_d1 = find(xi == -1);
fron_n = find(xi == 9);
fron_d = setdiff(fron_tot,fron_n);

%     Comprobación frontera
%     figure
%     hold on
%     scatter(xi(fron_tot),yi(fron_tot),'b')
%     scatter(xi(fron_d),yi(fron_d),'r','.')
%     hold off

% fron_d = setdiff(fron_tot,fron_n)'; %transponer aqui importante

fron_d = [fron_d fron_d+length(xi)];


%A = [nu*R 0*R -Cxl; 0*R nu*R -Cyl; Cxl' Cyl' zeros(size(Cxl,2))]; 
A = [nu*R 0*R -(1/rho)*Cxl; 0*R nu*R -(1/rho)*Cyl; Cxl' Cyl' zeros(size(Cxl,2))];

A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i) = 1;
end    

g1i = 0*xi';
g2i = 0*xi';
g1i(fron_d1) = 1-yi(fron_d1).^2; 
gg = [g1i;g2i;zeros(size(Cxl,2),1)];


f1i = f1(xi)';
f2i = f2(xi)';

vect_b = [M*f1i; M*f2i; zeros(size(Cxl,2),1)];
vect_b = vect_b - A*gg;

vect_b(fron_d)=0;
wh = A0\vect_b;
uh = wh + gg;

uh1 = uh(1:length(xi));
uh2 = uh((1:length(xi))+length(xi));
ph = uh((2*length(xi)+1):end);  
% Presion en nodos lineales, velocidad en cuadraticos. No son del mismo
% tamaño.

figure
trisurf(elem,xi,yi,uh1);
title('Vx')
figure
trisurf(elem,xi,yi,uh2);
title('Vy')
figure
quiver(xi,yi,uh1',uh2',0.7); % Genera campo de velocidades, el 0.7 escala
axis equal
figure
trisurf(elem,xi(1:size(Cxl,2)),yi(1:size(Cxl,2)),ph)
title('Presion')

% Especificamos los coeficientes de la PDE utilizando la forma general:
% m * d²u/dt² + d * du/dt - ∇·(c ∇u) + a * u = f
% - m = 0: No incluye términos de aceleración.
% - d = 0: No incluye términos de amortiguamiento.
% - c = 1: Coeficiente de difusión (gradiente en la ecuación).
% - a = 1: Multiplicador de la solución u.
% - f = 0: Sin término fuente.
% Estos valores definen una ecuación diferencial elíptica simple que se
% resuelve de manera estacionaria.
