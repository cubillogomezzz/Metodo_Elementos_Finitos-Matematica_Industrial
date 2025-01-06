% --- Ecuación del calor 1d parabólica discretizada en tiempo mediante
% euler explícito, condiciones dirichlet--- 


% Parámetros del problema, Ut + a*U - k*Uxx = f
f = @(x,t) 0*x + 0*t;
a=0;
k=0.1;

% Mallado temporal 
dt = 0.0001;
t0 = 0;
tf = 100;
Nt = (tf-t0)/dt;

% Mallado espacial
c=0;
d=1;
Ne = 100; % Número de elementos
h = (d-c)/Ne;
xi = c:h:d;

% Matrices de masas y rigideces, solo dependen de mallado espacial
ensamblaje_mat_masrig_1d_lin

% Matriz de rigidez global característica del problema
A = M; 

% Construcción A0
A0 = A;
A0([1 end],:)=0;
A0(:,[1 end])=0;
A0(1,1)=1;
A0(end,end)=1;

% Condiciones dirichlet 
gi = xi'*0; 
gi(1)=0;
gi(end)=0;

% Resolución iterativa

uhn = (sin(pi*xi')).^100; % Condición inicial
plot(xi,uhn);
title('0');
pause


for n = 1:Nt
    tn_mas1 = n*dt;
    fi = f(xi',(n-1)*dt); %f_n

    vect_b = M*uhn + dt*M*fi - dt*k*R*uhn - dt*a*M*uhn - A*gi; % Vector de carga característico del problema
    vect_b([1 end]) = 0; 

    whn = A0\vect_b;
    uhn = whn + gi; 

%     % Plot superpesto
%     hold on
%     plot(xi,uhn)
%     pause
%     hold off
    
    % Plot secuencial
    hold off 
    plot(xi,uhn,'b')
    title(sprintf('t = %.2f seconds', tn_mas1));
    pause(0.1)

end

hold off