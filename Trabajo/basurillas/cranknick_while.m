% Trabajo MEF Crank Nicolson

sistema = 0;

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

dg1 = @(x) at*(1-x/kt) + at*(-1/kt)*x;
dg2 = @(x) an*(1-x/kn) + an*(-1/kn)*x;

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

    T=0;
    N=0;

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

            ki = [al_tn*v(uhi(1:Nn),uhi(Nn+1:2*Nn)); al_nt*v(uhi(1:Nn),uhi(Nn+1:2*Nn))]; %bien

            % Calculo vector carga del tercio superior del sistema
           
            vec_carg = A*uhi - 0.5*dt*B*uhi + 0.5*dt*A*gi -0.5*dt*A*ki; %bien

            % Newton Raphson para sistemas de ecuaciones 
              
            % Inicialización

            Xn=[uhi;ki;gi]; 
            
            tol=1e-6;      % desviación relativa entre dos ultimas estimaciones
            N=400;        % el parámetro de finalización debería de ser tol, N grandes
            
            i=0;
            Xn_ant_prenew = Xn;
            while 1 % Newton
                i=i+1;
                Xn_ant = Xn;
                salida = 0;
                % Construcción F,JF

                F_Xn_ter1 = A*Xn(1:2*Nn) - 0.5*dt*A*Xn(4*Nn+1:6*Nn) + 0.5*dt*B*Xn(1:2*Nn) + 0.5*dt*A*Xn(2*Nn+1:4*Nn) - vec_carg;
                aux_g = sparse(2*Nn);
                F_Xn_ter2  = zeros(2*Nn,1);
                F_Xn_ter3  = zeros(2*Nn,1);

                di = eye(2*Nn);
                ze2 = sparse(2*Nn,2*Nn);
                eye_u1 = sparse(Nn,Nn);
                eye_u2 = sparse(Nn,Nn);% para construccion aux_k
                
                ze1 = sparse(Nn,Nn);

                for s=1:Nn %bucle auxiliar para construccion de F, DF

                    F_Xn_ter2(s) = Xn(2*Nn+s) - al_tn*Xn(s)*Xn(s+Nn); %bien
                    F_Xn_ter2(s+Nn) = Xn(2*Nn+s+Nn) - al_nt*Xn(s)*Xn(s+Nn);

                    F_Xn_ter3(s) = Xn(4*Nn+s) - g1(Xn(s)); %bien
                    F_Xn_ter3(s+Nn) = Xn(4*Nn+s+Nn) - g2(Xn(s+Nn));

                    aux_g(s,s) = -dg1(Xn(s)); %bien
                    aux_g(s+Nn,s+Nn) = -dg2(Xn(s+Nn));

                    eye_u1(s,s) = Xn(s); %bien
                    eye_u2(s,s) = Xn(s+Nn);
                end 

                F_Xn = [F_Xn_ter1; F_Xn_ter2; F_Xn_ter3];

%               aux_k = [(-1)*al_tn*eye_u2 ze1;ze1 (-1)*al_nt*di1] + [ze1 (-1)*al_tn*di1; (-1)*al_nt*di1 ze1];
                aux_k = [(-1)*al_tn*eye_u2 (-1)*al_tn*eye_u1; (-1)*al_nt*eye_u2 (-1)*al_nt*eye_u1]; %bien
                

                JF_Xn = [(A+0.5*dt*B) 0.5*dt*A 0.5*(-1)*dt*A;aux_k di ze2;aux_g ze2 di]; %bien
                
                paso = JF_Xn\F_Xn;
                normapaso = norm(paso)
                fprintf('pause newton')
                pause
                Xn=Xn-paso;
                dist = norm(Xn-Xn_ant);
                if dist/norm(Xn) < tol 
                    dist
                    salida = 1
                    i
                    break;
                end

                % actualizar cosas?
            
            end
            norm_actualiz = norm(Xn-Xn_ant_prenew)
            
            T=0;
            N=0;

            % Representacion gráfica

            if mod(n,10) == 0       
                     uhi = Xn(1:2*Nn);


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

                     
                     fprintf('pause tiempo')
                     pause
           end

end


% Analisis

paso(1:2*Nn)
paso(2*Nn+1:4*Nn)
paso(4*Nn+1:6*Nn)




