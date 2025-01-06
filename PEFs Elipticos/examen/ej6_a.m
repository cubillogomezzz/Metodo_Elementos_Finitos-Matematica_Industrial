% --- Integración numérica sobre un mallado de triángulos importado usando centro---

f = @(x,y) sin(y).*x.*exp(x-(1/9)*y.^2)/(x+y+1);

Q=0;
m = 1/3;
for i=1:length(t) 

% Indices de nodos que forman el elemento
n1=t(1,i);
n2=t(2,i);
n3=t(3,i);

% Transformación
A=[p(1,n2)-p(1,n1) p(1,n3)-p(1,n1);p(2,n2)-p(2,n1) p(2,n3)-p(2,n1)];
F = @(x,y) A*[x;y]+[p(1,n1);p(1,n1)];

w1=0.5;
p1 = F(m,m);
Q=Q+det(A)*w1*f(p1(1),p1(2));

end
Q