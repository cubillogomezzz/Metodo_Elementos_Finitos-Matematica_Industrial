% ---Integración numerica en mallado de cuadrilateros mediante distintas cuadraturas---

f=@(x,y) 1 + x*0;

% Funciones base bilinieales para transformación
phi1=@(x,y) 0.25*(y-1)*(x-1);
phi2=@(x,y) -0.25*(y-1)*(x+1);
phi3=@(x,y) 0.25*(y+1)*(x+1);
phi4=@(x,y) -0.25*(y+1)*(x-1);

Qt=0;
Q1p=0;
Qs=0;

for i=1:4

    % Coordenadas vértices del elemento 
    a1=[x(i,1);x(i,2)];  % modificar para formato pet
    a2=[x(i,3);x(i,4)];
    a3=[x(i,5);x(i,6)];
    a4=[x(i,7);x(i,8)];

    % Transformación y jacobiano
    F = @(x, y) [a1(1) * phi1(x, y) + a2(1) * phi2(x, y) + a3(1) * phi3(x, y) + a4(1) * phi4(x, y);
                  a1(2) * phi1(x, y) + a2(2) * phi2(x, y) + a3(2) * phi3(x, y) + a4(2) * phi4(x, y)];
    JF = @(x, y) [a1(1) * 0.25 * (y - 1) - a2(1) * 0.25 * (y - 1) + a3(1) * 0.25 * (y + 1) - a4(1) * 0.25 * (y + 1), ...
               a1(1) * 0.25 * (x - 1) - a2(1) * 0.25 * (x + 1) + a3(1) * 0.25 * (x + 1) - a4(1) * 0.25 * (x - 1); ...
               a1(2) * 0.25 * (y - 1) - a2(2) * 0.25 * (y - 1) + a3(2) * 0.25 * (y + 1) - a4(2) * 0.25 * (y + 1), ...
               a1(2) * 0.25 * (x - 1) - a2(2) * 0.25 * (x + 1) + a3(2) * 0.25 * (x + 1) - a4(2) * 0.25 * (x - 1)];
    dtm = @(x, y) det(JF(x, y)); % La transformación dilata de forma distinta las áreas en distitnas zonas del elemento

    % Coordenadas resto de puntos relevantes para la cuadratura 
    a5=F(0,-1); 
    a6=F(1,0);
    a7=F(0,1);
    a8=F(-1,0);
    a9=F(0,0);

    % Simpson
    Qs = Qs + (1/9) * (f(a1(1), a1(2)) * dtm(-1, -1) + f(a2(1), a2(2)) * dtm(-1, 1) + f(a3(1), a3(2)) * dtm(1, -1) + f(a4(1), a4(2)) * dtm(1, 1) ...
            + 4 * (f(a5(1), a5(2)) * dtm(0, -1) + f(a6(1), a6(2)) * dtm(1, 0) + f(a7(1), a7(2)) * dtm(0, 1) + f(a8(1), a8(2)) * dtm(-1, 0)) ...
            + 16 * f(a9(1), a9(2)) * dtm(0, 0));

    % 1 punto
    Q1p = Q1p + 4*f(a5(1),a5(2))*dtm(0,0);

    % Trapecios 
    Qt = Qt + f(a1(1), a1(2)) * dtm(-1, -1) + f(a2(1), a2(2)) * dtm(-1, 1) + f(a3(1), a3(2)) * dtm(1, -1) + f(a4(1), a4(2)) * dtm(1, 1);

end

Q1p
Qt
Qs