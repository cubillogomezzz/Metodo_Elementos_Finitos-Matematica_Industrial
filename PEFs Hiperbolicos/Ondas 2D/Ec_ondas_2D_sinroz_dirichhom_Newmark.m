% --- Problema hiperbólico ec ondas 2D con Newmark sin rozamiento ---
% Condiciones dirchlet homogéneas, beta=!0

% utt -c^2*uxx = f
% u(x,0) = h0(x)
% ut(x,0) = j0(x)


% Definición problema
c_cuad = 1;  
f = @(x,y,t) 0*x;

h0 = @(x,y) 0.5*exp(-((x-1).^2+(y-0.3).^2)/0.01);
j0 = @(x,y) x*0;

% Parámetros Newmark
gamma = 0.5; % 0<gamma<1 
beta = 0.25; % 0<beta<0.5, explícito si = 0

% Discretización espacial, mallado importado
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

fron_tot = unique([e(1,:) e(2,:)]);

% Discretización temporal 
t0=0;
tf=2.1;

Ne_t = 10000;
dt = (tf-t0)/Ne_t;

% Ensamblaje matrices
[R, M]=assema(p,t,1,1,0);

% Condiciones contorno para mat. rig. global de elíptico principal

A = M + beta*c_cuad*dt^2*R; % Caracterísitco de problema y esquema en t

A0 = A;
A0(fron_tot,:)=0;
A0(:,fron_tot)=0;
for i=fron_tot
    A0(i,i)=1;
end    

% Resolución 

    % Inicialización, obtención a en t=0
      fi = f(xi,yi,0)';
      ui = h0(xi,yi)';
      vi = j0(xi,yi)';  

      ai = M\(M*fi-c_cuad*R*ui); %(R sin imponerle condiciones, ya lo hace ui)

      figure(1)
      trisurf(elem,xi,yi,ui)
      title(0)
      axis([-1 1 -1 1 -1 1])
      pause

    % Iteración en tiempo
      for i=1:Ne_t

            % Elíptico principal

            ui_ant = ui;
            vi_ant = vi;
            ai_ant = ai;
            fi_ant = fi;

            fi = f(xi,yi,i*dt)';
           
            vec_carg = M*(ui_ant + dt*vi_ant + 0.5*dt^2*(1-2*beta)*ai_ant + beta*dt^2*fi); % Caracrterístico de problema y esquema
            vec_carg(fron_tot) = 0;

            ui = A0\vec_carg;

            % Obtención a, v en n+1 (en este caso no es necesario
            % martillazo mef)

            ai = (1/(beta*dt^2))*(ui-ui_ant-dt*vi_ant-0.5*dt^2*(1-2*beta)*ai_ant);
            vi = vi_ant + dt*(gamma*ai+(1-gamma)*ai_ant);

            % Representacion gráfica

                if mod(i,10) == 0
                     figure(1)
                     trisurf(elem,xi,yi,ui)
                     title(i*dt)
                     axis([-0.5 1.5 -0.5 1.5 -1 1])
                     %view(2)
                     shading interp
                     caxis([-.02 .02])
                     pause(0.001)
                end

      end



%---------- PROBLEMA HIPERBÓLICO -----------%
%
% utt - 2*k*ut - c2*lapl(u) = f
%
% u(gamma_d) = 0
% u(x,0) = u0
% ut(x,0) = v0
%
%------- discretización: NEWMARK --------%
%
% suponemos k=0, beta!=0
%
% a_0 = c2*lapl(u_0) + f_0
% u_n+1 = u_n + dt*v_n + 0.5*dt^2*(2*beta*a_n+1 + (1-2beta)*a_n)
% a_n+1 se despeja de lo anterior
% v_n+1 = v_n + dt*(gamm*a_n+1 + (1-gamma)*a_n)
%
%------- forma matricial --------%
%
% A = M  + (beta*c2*dt^2)R
%
% M*a_0 = -R*u_0 + M*f_0
% A0*u_n+1 = M(u_n + dt*v_n + 0.5*dt^2*(1-2*beta)*a_n + beta*dt^2*f_n+1)
% beta*dt^2 * a_n+1 = u_n+1 - u_n - dt*v_n - 0.5*dt^2*(1-2*beta)*a_n)
% v_n+1 = v_n + dt*(gamma*a_n+1 + (1-gamma)*a_n)



