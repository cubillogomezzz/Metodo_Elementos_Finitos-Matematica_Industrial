% --- Interpolación en R2 mediante elemento de referencia triangular con polinomio grado 2 --- 

% Para un elemento conocidos sus vertices y la función en estos 

f = @(x,y) sin(pi*x).*sin(pi*y) ;

X = [0.32; 0.74]; % Punto de estudio 


X1 = [0.25;0.5]; % Vértices del elemento
X2 = [0.5;0.75];
X3 = [0.25;0.75];
X4 = [0.375;0.625];
X5 = [0.375;0.75];
X6 = [0.25;0.625];

f1 = f(X1(1),X1(2)); 
f2 = f(X2(1),X2(2));
f3 = f(X3(1),X3(2));
f4 = f(X4(1),X4(2));
f5 = f(X5(1),X5(2));
f6 = f(X6(1),X6(2));

phi1= @(x,y) (1-x-y).*(1-2*x-2*y); % Funciones base de la interpolación
phi2 = @(x,y) x.*(2*x-1);
phi3 = @(x,y) y.*(2*y-1);
phi4 = @(x,y) (2-2*x-2*y).*(2*x);
phi5 = @(x,y) (2*y).*(2*x);
phi6 = @(x,y) (2-2*x-2*y).*(2*y);

A = [X2-X1 X3-X1]; % Transformacion afín desde elemento de referencia 
                   % (parte afín es + X1)

Xg = A\(X-X1); % Transformación inversa

xg = Xg(1);
yg = Xg(2);

interpol = f1*phi1(xg,yg) + f2*phi2(xg,yg) + f3*phi3(xg,yg) + f4*phi4(xg,yg) + f5*phi5(xg,yg) + f6*phi6(xg,yg) 
%(hacer más elegante con for si procede)

reales = f(X(1),X(2))
