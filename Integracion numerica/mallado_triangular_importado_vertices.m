% --- Integración numérica sobre un mallado de triángulos importado usando vértices---

f = @(x,y) 1+0*y;

Q=0;
for i=1:length(t) 

% Indices de nodos que forman el elemento
n1=t(1,i);
n2=t(2,i);
n3=t(3,i);

% Transformación
A=[p(1,n2)-p(1,n1) p(1,n3)-p(1,n1);p(2,n2)-p(2,n1) p(2,n3)-p(2,n1)];
F = @(x,y) A*[x;y]+[p(1,n1);p(1,n1)];

w1=0.5*(1/3);
w2=w1;
w3=w1;

p1 = F(0,0);
p2 = F(1,0);
p3 = F(0,1);

Q=Q+det(A)*(w1*f(p1(1),p1(2))+w2*f(p2(1),p2(2))+w3*f(p3(1),p3(2))); 

end
Q