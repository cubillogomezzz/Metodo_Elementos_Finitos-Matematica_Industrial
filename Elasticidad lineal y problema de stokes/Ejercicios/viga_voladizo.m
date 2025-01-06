% --- Ejercicio elasticidad lineal viga en voladizo ---

% -div(sigma(u))=fm
% constitutivas
% definicion tensor deformación
% condiciones


xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_tot = unique([e(1,:) e(2,:)]);
fron_d1 = find(xi == 11);
fron_d2 = find((xi > 10.5) & (yi == 0));
fron_d3 = find((xi > 10.5) & (abs(yi - 0.2) < 0.001));
fron_d = [fron_d1 fron_d2 fron_d3];


[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R21 M] = assema(p,t,[0 1 0 0]',1,0);
[R12 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);

E = 2.7e8;
mu = 0.2;
lambda = E*mu/((1+mu)*(1-2*mu));
nu = E/(2*(1+mu));
g = -9.8;
rho = 7900;

A11 = (lambda +2 *nu)*R11 + nu*R22;
A12 = lambda*R21 + nu*R12;
A21 = lambda*R12 + nu*R21;
A22 = (lambda +2 *nu)*R22 + nu*R11;

A = [A11 A12; A21 A22];

f1i = 0*xi';
f2i = 0*xi' + rho*g;

vect_b1 = M*f1i;
vect_b2 = M*f2i;

g1x = 0*xi'; %inicialización
g1y = 0*xi';

g2x = 0*xi';
g2y = 0*xi';

g3x = 0*xi';
g3y = 0*xi';

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
axis([-1 12 -0.5 0.7 0 1])

figure
trisurf(elem,xi+u1h',yi+u2h',0*xi)
axis([-1 12 -0.5 0.7 0 1])
view(2)



% Algunos detalles a tener en cuenta

% Tres tipos de condiciones: empotramiento, desplazamiento, fuerzas sup.
% nulas

% Al imponer desplazamientos puede quedar oscurecido el efecto de la
% rigidez

% Para calcular deformaciones, utilizar martillazo mef: hacer f =
% parcial(desplazamiento)/parcialx, truco del almendruco práctica 3