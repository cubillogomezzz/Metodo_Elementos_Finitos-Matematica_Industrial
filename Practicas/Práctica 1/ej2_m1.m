% --- Interpolación en R2 mediante elemento de referencia cuadrado con polinomio grado 2 --- 

% Para un elemento conocidos sus vertices y la función en estos 

f = @(x,y) sin(pi*x).*sin(pi*y) ;

[xg,yg]=meshgrid(-1:0.2:1,-1:0.2:1);% Punto de estudio 

X1 = [0 0.6 0.5 0;0 0 0.5 0.3]; % Nodos del elemento, entrada
X2 = [0.6 1 1 0.5;0 0 0.4 0.5];
X3 = [0.5 1 1 0.5;0.5 0.4 1 1];
X4 = [0 0.5 0.5 0;0.3 0.5 1 1];

for j=1:4

    nodes = [X1(:,j) X2(:,j) X3(:,j) X4(:,j)];
    fnode = f(nodes(1,:),nodes(2,:));
    
    % Funciones base para m=1, necesarias para transformación
    
    phi1_m1 = @(x,y) (0.5-0.5*x).*(0.5-0.5*y);
    phi2_m1 = @(x,y) (0.5+0.5*x).*(0.5-0.5*y);
    phi3_m1 = @(x,y) (0.5+0.5*x).*(0.5+0.5*y);
    phi4_m1 = @(x,y) (0.5-0.5*x).*(0.5+0.5*y);
    
    func_base = {phi1_m1, phi2_m1, phi3_m1, phi4_m1};
    
    % Transformación del elemento de referencia al resto de elementos
 

    Fx = @(x,y) phi1_m1(x,y).*nodes(1,1)+phi2_m1(x,y).*nodes(1,2)+... 
        phi3_m1(x,y).*nodes(1,3)+phi4_m1(x,y).*nodes(1,4);

    Fy = @(x,y) phi1_m1(x,y).*nodes(2,1)+phi2_m1(x,y).*nodes(2,2)+... 
        phi3_m1(x,y).*nodes(2,3)+phi4_m1(x,y).*nodes(2,4);

    xg_trans = Fx(xg,yg);
    yg_trans = Fy(xg,yg);

    % Evaluación de fh 

    interpol = 0*xg_trans;
    
    for i=1:4
    
        interpol = interpol + fnode(i).*func_base{i}(xg,yg);
     
    end

    hold on 
    surf(xg_trans,yg_trans,interpol)

end    
    

