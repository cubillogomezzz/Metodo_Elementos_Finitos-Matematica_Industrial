% --- Interpolación en R2 mediante elemento de referencia cuadrado con polinomio grado 2 --- 

% Para un elemento conocidos sus vertices y la función en estos 

f = @(x,y) sin(pi*x).*sin(pi*y) ;

X = [0.6; 0.24]; % Punto de estudio 


X1 = [0.25;0.5]; % Nodos del elemento, entrada
X2 = [0.5;0.5];
X3 = [0.5;0.75];
X4 = [0.25;0.75];
X5 = [0.625;0.5];
X6 = [0.5;0.625];
X7 = [0.375;0.75];
X8 = [0.25;0.625];
X9 = [0.375;0.625];

nodes = [X1 X2 X3 X4 X5 X6 X7 X8 X9];
fnode = f(nodes(1,:),nodes(2,:));

phi1 = @(x,y) (0.5-0.5*x).*(-x).*(-y).*(0.5-0.5*y); % Funciones base de la interpolación
phi2 = @(x,y) (0.5+0.5*x).*(x).*(-y).*(0.5-0.5*y);
phi3 = @(x,y) (0.5+0.5*x).*(x).*(y).*(0.5+0.5*y);
phi4 = @(x,y) (0.5-0.5*x).*(-x).*(y).*(0.5+0.5*y);

phi5 = @(x,y) (1-x).*(x+1).*(-y).*(0.5-0.5*y);
phi6 = @(x,y) (1-y).*(y+1).*(x).*(0.5+0.5*x);
phi7 = @(x,y) (1-x).*(x+1).*(y).*(0.5+0.5*y);
phi8 = @(x,y) (1-y).*(y+1).*(-x).*(0.5-0.5*x);

phi9 = @(x,y) (1-y).*(y+1).*(1-x).*(x+1);

func_base = {phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9};

% Funciones base para m=1, necesarias para transformación

phi1_m1 = @(x,y) (0.5-0.5*x).*(0.5-0.5*y);
phi2_m1 = @(x,y) (0.5+0.5*x).*(0.5-0.5*y);
phi3_m1 = @(x,y) (0.5+0.5*x).*(0.5+0.5*y);
phi4_m1 = @(x,y) (0.5-0.5*x).*(0.5+0.5*y);

% Funciones auxiliares para construir DF en newton

aux = @(t) (0.5+0.5*t)*0.5;

% Transformación al elemento de referencia 

F = @(x) [phi1_m1(x(1),x(2)).*X1(1)+phi2_m1(x(1),x(2)).*X2(1)+... % Transformación en forma de sistema no lineal 2x2;
    phi3_m1(x(1),x(2)).*X3(1)+phi4_m1(x(1),x(2)).*X4(1)-X(1); ...
    phi1_m1(x(1),x(2)).*X1(2)+phi2_m1(x(1),x(2)).*X2(2)+... 
    phi3_m1(x(1),x(2)).*X3(2)+phi4_m1(x(1),x(2)).*X4(2)-X(2)]; 

JF = @(x) [(-1)*aux(-x(2)).*X1(1)+aux(-x(2)).*X2(1)+...
    aux(x(2)).*X3(1)+aux(x(2))*(-1).*X4(1) (-1)*aux(-x(1)).*X1(1)+(-1)*aux(x(1)).*X2(1)+...
    aux(x(1)).*X3(1)+aux(-x(1)).*X4(1) ; 
    (-1)*aux(-x(2)).*X1(2)+aux(-x(2)).*X2(2)+...
    aux(x(2)).*X3(2)+aux(x(2))*(-1).*X4(2) (-1)*aux(-x(1)).*X1(2)+(-1)*aux(x(1)).*X2(2)+...
    aux(x(1)).*X3(2)+aux(-x(1)).*X4(2)]; % Jacobiana, cuidado espacios y ... 



Xn=[0;0]; % Conocido el intervalo de la solución, tomamos (0,0) como punto inicial por estar centrado en este. 

tol=1e-10;      
N=1000;        

for i=1:N
    Xn_ant = Xn;
    Xn=Xn-JF(Xn)\F(Xn); 
    if norm(Xn-Xn_ant)/norm(Xn) < tol 
        break;
    end
end

interpol = 0;

for i=1:9

    interpol = interpol + fnode(i).*func_base{i}(Xn(1),Xn(2));
 
end

interpol
reales = f(X(1),X(2))
error = abs(interpol-reales)
