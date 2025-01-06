% --- Ecuación del calor con advección 2d parabólica discretizada en tiempo mediante
% crank-nicolson, condiciones dirichlet y neuman homogeneas --- 

% CN es implícito, segundo orden en tiempo, incondicionalmente estable,
% pero puede presentar oscilaciones si no se da la relación adecuada entre
% k,dt,h. Típicamente k*dt<0.5*h^2.

% Parámetros del problema, Ut + a*U + b.grad(u) - k*lap(U) = f
f = @(x,y,t) 0*x + 0*y + 0*t;
a=0; % Termino de reaccion. a<0 -> exotermica
k=0.1;
b1 = 2;
b2 = 1;

% Mallado temporal 
dt = 0.1;
t0 = 0;
tf = 100;
Nt = (tf-t0)/dt;

% Mallado espacial
xi = p(1,:); 
yi = p(2,:);
elem = t(1:3,:)'; 

% Fronteras y condiciones
fron_tot = unique(e(1:2,:))';
tol = 0.001;
aux = find((yi + xi - 0.5) < tol );
fron_n = find(yi == 1 | xi == 1); %cuidado extremos de los intervalos
fron_d2 = setdiff(fron_tot,[aux fron_n]); 
fron_d1 = setdiff(fron_tot,[fron_d2 fron_n]);
fron_d = [fron_d1 fron_d2];

gi = 0*xi'; 
gi(fron_d1)=1; % Establecer valores en frontera según condiciones 
gi(fron_d2)=0;

%     Comprobación frontera
    figure
    hold on
    scatter(xi(fron_tot),yi(fron_tot),'b')
    scatter(xi(fron_d),yi(fron_d),'r','.')
    hold off

% Matrices de masas y rigideces, solo dependen de mallado espacial
matrices_rig_mas_conv_2d

% Matriz de rigidez global característica del problema
A = (1+a*dt*0.5)*M + 0.5*k*dt*R + dt*0.5*b1*Cx + dt*0.5*b2*Cy; 
% Construcción A0, matando fila_i+col_j a cero con diagonal a 1 para frontera dirichlet
A0 = A;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end    

% Resolución iterativa

uhn = 0*xi';
figure
trisurf(elem,xi,yi,uhn)
title('0');
pause

for n = 1:Nt
    tn_mas1 = n*dt;

    fi = f(xi',yi',(n-1)*dt);
    fi_m1 = f(xi',yi',n*dt);  

    vect_b = 0.5*dt*M*(fi+fi_m1) + (1-dt*0.5*a)*M*uhn - 0.5*dt*k*R*uhn - b1*0.5*dt*Cx*uhn - b2*0.5*dt*Cy*uhn ...
     - A*gi ; % Vector de carga característico del problema %Antes teniamos otra vaina con 1-a

    vect_b(fron_tot) = 0;

    whn = A0\vect_b;  
    uhn = whn + gi; 

    
    % Plot secuencial
    hold off 
    trisurf(elem,xi,yi,uhn)
    axis([0 1 0 1 0 2])
    title(sprintf('t = %.2f seconds', tn_mas1));
    pause(0.05)

end

hold off

% Tanto dirichlet (homogeneas o no, al cambiar a w estas últimas) como
% neuman homogeneas anulan el termino integrado en frontera del laplaciano.
% Las dirchlet imponen cero en el sistema, las neuman dejan libre. 