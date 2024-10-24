% --- Interpolación en R2 mediante elemento de referencia triangular con polinomio grado 1 --- 

% Para un elemento conocidos sus vertices y la función en estos 

f = @(x,y) sin(pi*x).*sin(pi*y) ;

X = [0.32 0.30; 0.74 0.79]; % Puntos de estudio 

X1 = [0.25;0.5]; % Vértices del elemento
X2 = [0.5;0.75];
X3 = [0.25;0.75];

f1 = f(X1(1),X1(2)); 
f2 = f(X2(1),X2(2));
f3 = f(X3(1),X3(2));

phi1= @(x,y) 1-x-y; % Funciones base de la interpolación
phi2 = @(x,y) x;
phi3 = @(x,y) y;

A = [X2-X1 X3-X1]; % Transformacion afín desde elemento de referencia 
                   % (parte afín es + X1)

Xg = A\(X-X1); % Transformación inversa

xg = Xg(:,1);
yg = Xg(:,2);

f1*phi1(xg,yg) + f2*phi2(xg,yg) + f3*phi3(xg,yg) 
