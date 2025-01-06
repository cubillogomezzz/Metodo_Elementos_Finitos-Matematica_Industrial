%  --- Elastodinámica con Newmark implícito, con rozamiento --- 

% rho*utt - div(sigma(u)) = f
% CF Neuman homo, dirichlet
% u(x,0) = h0(x)
% ut(x,0) = j0(x)

% Parámetros elasticidad
E = 2.7e8;
mu = 0.2;
lambda = E*mu/((1+mu)*(1-2*mu));
nu = E/(2*(1+mu));
g = -9.8;
rho = 7900;
umax = 0;
b = 20e5; %rozamiento

% Definición problema
f1 = @(x,y,t) 0*x;
f2 = @(x,y,t) 0*x+rho*g;

h01 = @(x,y) x*0;
h02 = @(x,y) x*0;
j01 = @(x,y) x*0;
j02 = @(x,y) x*0;

% Parámetros Newmark
gamma = 0.5; % 0<gamma<1 
beta = 0.25; % 0<beta<0.5, explícito si = 0

% Discretización espacial, mallado importado
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

% Fronteras y condiciones
fron_tot = unique([e(1,:) e(2,:)]);
fron_d1 = find(xi == 11);
fron_d2 = find((xi > 10.5) & (yi == 0));
fron_d3 = find((xi > 10.5) & (abs(yi - 0.2) < 0.001));
fron_d = [fron_d1 fron_d2 fron_d3];

%     %Comprobación frontera
%     figure
%     hold on
%     scatter(xi(fron_tot),yi(fron_tot),'b')
%     scatter(xi(fron_d),yi(fron_d),'r','.')
%     hold off

% g1x = 0*xi'; % Inicialización
% g1y = 0*xi';
% 
% g2x = 0*xi';
% g2y = 0*xi';
% 
% g3x = 0*xi';
% g3y = 0*xi';

% g1x(fron_d1) = (1-cos(pi/4))*(1-xi(fron_d1)'); % Establecer condiciones
% g1y(fron_d1) = sin(pi/4)*(1-xi(fron_d1)'); 

% g2x(fron_d2) = (1-cos(pi/4))*(1-xi(fron_d2)');
% g2y(fron_d2) = (-1)*sin(pi/4)*(1-xi(fron_d2)');

% Discretización temporal 
t0=0;
tf=30;

Ne_t = 100000;
dt = (tf-t0)/Ne_t;

% Ensamblaje matrices
[R11 M] = assema(p,t,[1 0 0 0]',1,0);
[R21 M] = assema(p,t,[0 1 0 0]',1,0);
[R12 M] = assema(p,t,[0 0 1 0]',1,0);
[R22 M] = assema(p,t,[0 0 0 1]',1,0);

A11 = (lambda +2 *nu)*R11 + nu*R22;
A12 = lambda*R21 + nu*R12;
A21 = lambda*R12 + nu*R21;
A22 = (lambda +2 *nu)*R22 + nu*R11;

A = [A11 A12; A21 A22];
aux_cero = sparse(size(M,1),size(M,1));
M_barra = [M aux_cero;aux_cero M];

% Sistema característico 
Q = 1 + dt*gamma*b*0.5*rho;
W = 1-(dt*gamma*b*0.5)/(Q*rho);
Wb = 1-(dt*beta*b*0.5)/(Q*rho);

Aglob = M_barra + beta*dt.^2*(1/rho)*W*A ; 

% Condiciones contorno para mat. rig. global de elíptico principal
fron_d = [fron_d fron_d+length(xi)]; % Al aplanar vector, doblar frontera

A0 = Aglob;
A0(fron_d,:)=0;
A0(:,fron_d)=0;
for i=fron_d
    A0(i,i)=1;
end    

% Resolución 

    % Inicialización, obtención a en t=0

      fi1 = f1(xi,yi,0)';
      fi2 = f2(xi,yi,0)';
      fi = [fi1;fi2];

      ui1 = h01(xi,yi)';
      ui2 = h02(xi,yi)';
      ui = [ui1;ui2];

      vi1 = j01(xi,yi)';
      vi2 = j02(xi,yi)';
      vi = [vi1;vi2]; 

      ai = M_barra\((M_barra*(fi-0.5*b*vi)-A*ui)*(1/rho)); 

      figure(1)
      trisurf(elem,xi+ui1',yi+ui2',0*xi)
      title(0)
      view(2)
      axis([-1 12 -0.5 0.7 0 1])
      pause

   
    % Iteración en tiempo
      for i=1:Ne_t

            % Elíptico principal

            ui_ant = ui;
            vi_ant = vi;
            ai_ant = ai;
            fi_ant = fi;

            fi1 = f1(xi,yi,dt*i)';
            fi2 = f2(xi,yi,dt*i)';
            fi = [fi1;fi2];
           

            vec_carg = M_barra*(ui_ant + dt*Wb*vi_ant + dt^2*beta*(1/rho)*W*fi + ...
                ((1-2*beta)*0.5*dt.^2-(1-gamma)*beta*dt.^2*(0.5*dt*b/(Q*rho)))*ai_ant);
            % Caracrterístico de problema y esquema
            vec_carg(fron_d) = 0;

            ui = A0\vec_carg;
            
            % Obtención a, v en n+1 

            ai = 1/(beta*dt^2)*(ui-ui_ant-dt*vi_ant-0.5*dt^2*(1-2*beta)*ai_ant);
            vi = vi_ant + dt*(gamma*ai+(1-gamma)*ai_ant);

            ext = find(xi == 0 & yi == 0);
            u_ext = abs(ui(ext+length(xi)));
            if u_ext > umax
                umax=u_ext;
            end    
            % Representacion gráfica

                if mod(i,10) == 0              
                    u1h = ui(1:length(xi));
                    u2h = ui(length(xi)+1 : end);
                    trisurf(elem,xi+u1h',yi+u2h',0*xi)
                    axis([-1 12 -0.5 0.7 0 1])
                    title(i*dt)
                    view(2)
                    pause(0.01)
                end

      end



%------- discretización: NEWMARK --------%
%
% suponemos k=0, beta!=0
%
% a_0 = c2*lapl(u_0) + f_0
% u_n+1 = u_n + dt*v_n + 0.5*dt^2*(2*beta*a_n+1 + (1-2beta)*a_n)
% a_n+1 se despeja de lo anterior
% v_n+1 = v_n + dt*(gamm*a_n+1 + (1-gamma)*a_n)