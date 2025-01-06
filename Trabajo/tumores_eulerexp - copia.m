% Trabajo MEF Crank Nicolson

an = 5;   %[1/día]
at = 0;    %[1/día]
Dn = 2.0e-15;   %[cm^2/día]
Dt = 4.2e-3;    %[cm^2/día]
kn = 1; %capacidad carga del medio
kt = 1;
al_tn = 1.0e-4;   %[1/día]
al_nt = 0.2;   %[1/día

h0 = @(x,y) exp(-((x-0.5).^2+(y-0.5).^2)/0.01); %condicion inicial

g1 = @(x) at*(1-x/kt).*x; % funciones término reacción
g2 = @(x) an*(1-x/kn).*x; 

v = @(x,y) x.*y; % función términos competición

t0 = 0;
tf = 10;
Ne_t = 1000;
dt = (tf-t0)/Ne_t;


xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

Nn = length(xi);

fron_tot = unique(e(1:2,:));

[R M]=assema(p,t,1,1,0);

A = [M 0*M; M*0 M]; % mat rig. global, sin cond dirich
B = [Dt*R 0*R; 0*R Dn*R];


T0 = h0(xi,yi)';
N0 = (0*T0+1);
uhi = [T0; N0]; %uhi_0

% Integral para cantidad total
% 
%     T=0;
%     N=0;
%     for i=1:length(t) 
%     
%     n1=t(1,i);
%     n2=t(2,i);
%     n3=t(3,i);
% 
%     Atra=[p(1,n2)-p(1,n1) p(1,n3)-p(1,n1);p(2,n2)-p(2,n1) p(2,n3)-p(2,n1)];
%     F = @(x,y) Atra*[x;y]+[p(1,n1);p(1,n1)];
%     
%     w1=0.5*(1/3);
%     w2=w1;
%     w3=w1;
%     
%     p1 = F(0,0);
%     p2 = F(1,0);
%     p3 = F(0,1);
%     
%     T=T+det(Atra)*(w1*uhi(n1)+w2*uhi(n2)+w3*uhi(n3)); 
%     N=N+det(Atra)*(w1*uhi(Nn + n1)+w2*uhi(Nn + n2)+w3*uhi(Nn + n3)); 
%     
%     end
    
    
% Plot t=0

figure(1)

subplot(2, 1, 1)
trisurf(elem,xi,yi,uhi(1:Nn))
%view(2)
title(['cantidad total T = ', num2str(T)]);
axis([0 1 0 1 0 1])
colorbar
caxis([0 1])

subplot(2, 1, 2)
trisurf(elem,xi,yi,uhi(Nn+1:2*Nn))
%view(2)
title(['cantidad total T = ', num2str(N)]);
axis([0 1 0 1 0 1])
colorbar
caxis([0 1])

sgtitle(['tiempo = ', num2str(0)]);

pause

for n=1:Ne_t

            % Calculo g, k iteración n

            gi1 = g1(uhi(1:Nn));
            gi2 = g2(uhi((Nn+1):2*Nn));
            gi = [gi1; gi2];

            ki = [al_t*v(uhi(1:Nn),uhi(Nn+1:2*Nn)); al_n*v(uhi(1:Nn),uhi(Nn+1:2*Nn))];

            % Calculo vector carga del tercio superior del sistema
           
            vec_carg = A*uhi - 0.5*dt*B*uhi + 0.5*dt*A*gi -0.5*dt*A*ki; 

            % Newton Raphson para sistemas de ecuaciones 
              
            % Inicialización

            Xn=[uhi;ki;gi]; 
            
            tol=1e-8;      % desviación relativa entre dos ultimas estimaciones
            N=1000;        % el parámetro de finalización debería de ser tol, N grandes
            
            for i=1:N
                Xn_ant = Xn;

                % Construcción F,JF

                F_Xn_ter1 = A*Xn(1:2*Nn) - 0.5*dt*A*Xn(4*Nn+1:6*Nn) + 0.5*dt*B*Xn(1:2*Nn) + 0.5*dt*A*Xn(2*Nn+1:4*Nn) - vec_carg;
                aux_g = sparse(2*Nn);

                for p=1:Nn
                    F_Xn_ter2  = zeros(2*Nn)';
                    F_Xn_ter2(p) = Xn(2*Nn+p) - al_tn*Xn(p)*Xn(p+Nn);
                    F_Xn_ter2(p+Nn) = Xn(2*Nn+p+Nn) - al_nt*Xn(p)*Xn(p+Nn);

                    F_Xn_ter3  = zeros(2*Nn)';
                    F_Xn_ter3(p) = Xn(4*Nn+p) - g1(Xn(p));
                    F_Xn_ter3(p+Nn) = Xn(4*Nn+p+Nn) - g2(Xn(p+Nn));

                    aux_g(p,p) = -g1(Xn(p));
                    aux_g(p+Nn,p+Nn) = -g2(Xn(p+Nn));
                end 

                F_Xn = [F_Xn_ter1; F_Xn_ter2; F_Xn_ter3];

                di = diag(2*Nn);
                ze = sparse(2*Nn);
                di1 = diag(Nn);
                ze1 = sparse(Nn)
                aux_k = [-al_tn*di1 ze1;ze1 -al_nt*di1] + [ze1 -al_tn*di1;-al_nt*di1 ze1];
                

                JF_Xn = [A+0.5*dt*B 0.5*dt*A -0.5*dt*A;aux_k diag ze;aux_g diag ze];

                Xn=Xn-JF(Xn)\F(Xn);
                if norm(Xn-Xn_ant)/norm(Xn) < tol 
                    break;
                end
            
               % hold on 
               % plot([Xn_ant(1) Xn(1)],[Xn_ant(2) Xn(2)], 'r.-') 
               % pause  %dibuja recta entre iteraciones, hacer sobre plot inicial
            
            end

            




            T=0;
            N=0;



            % Representacion gráfica

                if mod(n,10) == 0

%                     for j=1:length(t) 
%                     
%                     n1=t(1,j);
%                     n2=t(2,j);
%                     n3=t(3,j);
%                 
%                     Atra=[p(1,n2)-p(1,n1) p(1,n3)-p(1,n1);p(2,n2)-p(2,n1) p(2,n3)-p(2,n1)];
%                     F = @(x,y) Atra*[x;y]+[p(1,n1);p(1,n1)];
%                     
%                     w1=0.5*(1/3);
%                     w2=w1;
%                     w3=w1;
%                     
%                     p1 = F(0,0);
%                     p2 = F(1,0);
%                     p3 = F(0,1);
%                     
%                     T=T+det(Atra)*(w1*uhi(n1)+w2*uhi(n2)+w3*uhi(n3)); 
%                     N=N+det(Atra)*(w1*uhi(Nn + n1)+w2*uhi(Nn + n2)+w3*uhi(Nn + n3)); 
%                     
%                     end                    

                     figure(2)

                     subplot(2, 1, 1)
                     trisurf(elem,xi,yi,uhi(1:Nn))
                     axis([0 1 0 1 0 1])
                     title(['cantidad total T = ', num2str(T)]);
                     %view(2)
                     shading interp
                     colorbar
                     caxis([0 1])

                     subplot(2, 1, 2)
                     trisurf(elem,xi,yi,uhi(Nn+1:2*Nn))
                     axis([0 1 0 1 0 1])
                     title(['cantidad total N = ', num2str(N)]);
                     %view(2)
                     shading interp
                     colorbar
                     caxis([0 1])

                     sgtitle(['tiempo = ', num2str(dt*n)]);

                     pause(0.001)
                end

end
