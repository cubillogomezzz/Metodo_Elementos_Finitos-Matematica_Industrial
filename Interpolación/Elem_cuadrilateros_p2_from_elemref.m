% --- Interpolación en R2 mediante elemento de referencia cuadrado con polinomio grado 2 --- 
% Coordenadas dato en el elemento de referencia

% Para un elemento conocidos sus vertices y la función en estos 

f = @(x,y) sin(pi*x).*sin(pi*y) ;

% Coordenadas dato en elemento de referencia

%[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1); % Mallado
xg = ; % Un solo punto
yg = ;

X1 = [0 0.6 0.5 0;0 0 0.5 0.3]; % Nodos del elemento, entrada
X2 = [0.6 1 1 0.5;0 0 0.4 0.5];
X3 = [0.5 1 1 0.5;0.5 0.4 1 1];
X4 = [0 0.5 0.5 0;0.3 0.5 1 1];

X5 = 0.5 * (X1+X2);
X6 = 0.5 * (X2+X3);
X7 = 0.5 * (X3+X4);
X8 = 0.5 * (X4+X1);

X9 = 0.25*(X1+X2+X3+X4);

for j=1:4

    nodes = [X1(:,j) X2(:,j) X3(:,j) X4(:,j) X5(:,j) X6(:,j) X7(:,j) X8(:,j) X9(:,j)];
    fnode = f(nodes(1,:),nodes(2,:)); 
    
    % Funciones base para m=1, necesarias para transformación
    
    phi1_m1 = @(x,y) (0.5-0.5*x).*(0.5-0.5*y);
    phi2_m1 = @(x,y) (0.5+0.5*x).*(0.5-0.5*y);
    phi3_m1 = @(x,y) (0.5+0.5*x).*(0.5+0.5*y);
    phi4_m1 = @(x,y) (0.5-0.5*x).*(0.5+0.5*y);

    % Funciones base para m=2

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
    
    % Transformación del elemento de referencia al resto de elementos
 

    Fx = @(x,y) phi1_m1(x,y).*nodes(1,1)+phi2_m1(x,y).*nodes(1,2)+... 
        phi3_m1(x,y).*nodes(1,3)+phi4_m1(x,y).*nodes(1,4);

    Fy = @(x,y) phi1_m1(x,y).*nodes(2,1)+phi2_m1(x,y).*nodes(2,2)+... 
        phi3_m1(x,y).*nodes(2,3)+phi4_m1(x,y).*nodes(2,4);

    xg_trans = Fx(xg,yg);
    yg_trans = Fy(xg,yg);

    % Evaluación de fh 

    interpol = 0*xg_trans;
    
    for i=1:9
    
        interpol = interpol + fnode(i).*func_base{i}(xg,yg);
     
    end

    hold on 
    surf(xg_trans,yg_trans,interpol)

end    