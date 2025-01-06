f=@(X) sin(pi*X(1)).*sin(pi*X(2));
x=zeros(4,8);
%nodos de cada elemento
x(1,:)=[  0   0   0.6 0   0.5 0.5      0 0.3]; %coordenadas de los cuatro nodos del elemnto 1
x(2,:)=[0.6   0     1 0     1 0.4    0.5 0.5];
x(3,:)=[  1 0.4     1 1   0.5   1    0.5 0.5];
x(4,:)=[0.5   1     0 1     0 0.3    0.5 0.5];

Qs=0;
for i=1:4
    a1=[x(i,1);x(i,2)];
    a2=[x(i,3);x(i,4)];
    a3=[x(i,5);x(i,6)];
    a4=[x(i,7);x(i,8)];

    %Tengo que sacar las coordenadas de los nodos 5,6 ,7 ,8 y 9 del elemento de
    %referencia en elemento real
    %uso Fi para ello
    phi1=@(x,y) 0.25*(y-1)*(x-1);
    phi2=@(x,y) -0.25*(y-1)*(x+1);
    phi3=@(x,y) 0.25*(y+1)*(x+1);
    phi4=@(x,y) -0.25*(y+1)*(x-1);

    F=@(X) [x(i,1)*phi1(X(1),X(2))+ x(i,3)*phi2(X(1),X(2))+ x(i,5)*phi3(X(1),X(2))+ x(i,7)*phi4(X(1),X(2));
        x(i,2)*phi1(X(1),X(2))+ x(i,4)*phi2(X(1),X(2))+ x(i,6)*phi3(X(1),X(2))+ x(i,8)*phi4(X(1),X(2))];
    %los ptos en los que se evalúa F están en el elemnto de referencia. Las a
    %en el elemnento real
    a5=F([0,-1]); %el input de F es un vector
    a6=F([1,0]);
    a7=F([0,1]);
    a8=F([-1,0]);
    a9=F([0,0]);

    %la X es X_gorro
    JF=@(X) [x(i,1)*0.25*(X(2)-1)-x(i,3)*0.25*(X(2)-1)+x(i,5)*0.25*(X(2)+1)-x(i,7)*0.25*(X(2)+1)    x(i,1)*0.25*(X(1)-1)- x(i,3)*0.25*(X(1)+1)+ x(i,5)*0.25*(X(1)+1)- x(i,7)*0.25*(X(1)-1);
        x(i,2)*0.25*(X(2)-1)-x(i,4)*0.25*(X(2)-1)+x(i,6)*0.25*(X(2)+1)-x(i,8)*0.25*(X(2)+1)    x(i,2)*0.25*(X(1)-1)- x(i,4)*0.25*(X(1)+1)+ x(i,6)*0.25*(X(1)+1)- x(i,8)*0.25*(X(1)-1)];
    dtm=@(X) det(JF(X));

    %la función q se integra es f_gorro(x_g)*det(x_g) ahora el determinante no
    %es cte no sale fuera tb se integra
    Qs=Qs+(1/9)*(f(a1)*dtm([-1,-1])+f(a2)*dtm([-1,1])+f(a3)*dtm([1,-1])+f(a4)*dtm([1,1])+4*(f(a5)*dtm([0,-1])+f(a6)*dtm([1,0])+f(a7)*dtm([0,1])+f(a8)*dtm([-1,0]))+16*f(a9)*dtm([0,0]));
end
Qs
valor_real=4*pi^(-2);
error=abs(valor_real-Qs);