%  --- Problema elasticidad lineal cond dirchlet generales, neuman
%  homogeneas ---

% -div(sigma(u))=fm
% constitutivas
% definicion tensor deformación
% condiciones frontera: dirichlet (desplazamiento impuesto), neuman homog
% (fuerza superficial nula)


xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_tot = unique(e(1:2,:));
fron_d1 = find(yi == 1);
fron_d2 = find(yi == -1);
fron_d3 = find((xi == 1));
fron_d = [fron_d1 fron_d2 fron_d3];

%     Comprobación frontera
%     figure
%     hold on
%     scatter(xi(fron_tot),yi(fron_tot),'b')
%     scatter(xi(fron_d),yi(fron_d),'r','.')
%     hold off


[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R21 M] = assema(p,t,[0 1 0 0]',1,0);
[R12 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);

E = 2.7e8;
nu = 0.2; % dato nu, se calcula mu, al reves

lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

g = 0;
rho = 1000;

A11 = (lambda +2 *mu)*R11 + mu*R22;
A12 = lambda*R21 + mu*R12;
A21 = lambda*R12 + mu*R21;
A22 = (lambda +2 *mu)*R22 + mu*R11;

A = [A11 A12; A21 A22];

f1i = 0*xi'; %fuente
f2i = 0*xi';

vect_b1 = M*f1i;
vect_b2 = M*f2i;

g1x = 0*xi'; %inicialización
g1y = 0*xi';

g2x = 0*xi';
g2y = 0*xi';

g3x = 0*xi';
g3y = 0*xi';

% Establecer condiciones
g1x(fron_d1) = (1-cos(pi/4))*(1-xi(fron_d1)');
g1y(fron_d1) = sin(pi/4)*(1-xi(fron_d1)'); 

g2x(fron_d2) = (1-cos(pi/4))*(1-xi(fron_d2)');
g2y(fron_d2) = (-1)*sin(pi/4)*(1-xi(fron_d2)');

% Comprobar condiciones
figure
quiver(xi(fron_d1)',yi(fron_d1)',g1x(fron_d1),g1y(fron_d1))
figure
quiver(xi(fron_d2)',yi(fron_d2)',g2x(fron_d2),g2y(fron_d2))

gh = [g1x; g1y] + [g2x; g2y];

vect_b = [vect_b1; vect_b2] - A * (gh);

fron_d = [fron_d fron_d+length(xi)]; %cambio al doblarse u en u1 | u2

A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i) = 1;
end    

vect_b(fron_d)=0;

wh = A0\vect_b;
uh = wh + gh;

u1h = uh(1:length(xi));
u2h = uh(length(xi)+1 : end);

hold off
% figure
% trisurf(elem,xi,yi,u1h)
% figure
% trisurf(elem,xi,yi,u2h)
%dibujito deformacion, dibujar primero malla y luego suma
figure

trisurf(elem,xi,yi,0*xi)
view(2)
axis([-1.2 1.2 -1.2 1.2 0 1])

figure
trisurf(elem,xi+u1h',yi+u2h',0*xi)
view(2)



% Algunos detalles a tener en cuenta

% Tres tipos de condiciones: empotramiento, desplazamiento, fuerzas sup.
% nulas

% Al imponer desplazamientos puede quedar oscurecido el efecto de la
% rigidez

% Para calcular deformaciones, utilizar martillazo mef: hacer f =
% parcial(desplazamiento)/parcialx, truco del almendruco práctica 3

