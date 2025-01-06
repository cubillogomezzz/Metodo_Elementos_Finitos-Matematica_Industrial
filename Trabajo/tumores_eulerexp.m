% Trabajo MEF euler exp

an = 5;   %[1/día]
at = 0;    %[1/día]
Dn = 2.0e-15;   %[cm^2/día]
Dt = 4.2e-3;    %[cm^2/día]
kn = 1; %capacidad carga del medio
kt = 1;
al_tn = 1.0e-4;   %[1/día]
al_nt = 0.2;   %[1/día

h0 = @(x,y) exp(-((x-0.5).^2+(y-0.5).^2)/0.01); %condicion inicial

g1 = @(x) (1-x/kt).*x; % funciones término reacción
g2 = @(x) (1-x/kn).*x; 

k = @(x,y) x.*y; % función términos competición

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
C = [at*M 0*M; 0*M an*M];
D = [al_tn*M 0*M; 0*M al_nt*M];

T0 = h0(xi,yi)';
N0 = (0*T0+1);
uhi = [T0; N0]; %uhi_0

% Integral para cantidad total

    T=0;
    N=0;
    for i=1:length(t) 
    
    n1=t(1,i);
    n2=t(2,i);
    n3=t(3,i);

    Atra=[p(1,n2)-p(1,n1) p(1,n3)-p(1,n1);p(2,n2)-p(2,n1) p(2,n3)-p(2,n1)];
    F = @(x,y) Atra*[x;y]+[p(1,n1);p(1,n1)];
    
    w1=0.5*(1/3);
    w2=w1;
    w3=w1;
    
    p1 = F(0,0);
    p2 = F(1,0);
    p3 = F(0,1);
    
    T=T+det(Atra)*(w1*uhi(n1)+w2*uhi(n2)+w3*uhi(n3)); 
    N=N+det(Atra)*(w1*uhi(Nn + n1)+w2*uhi(Nn + n2)+w3*uhi(Nn + n3)); 
    
    end
    
    
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

for i=1:Ne_t

            
            gi1 = g1(uhi(1:Nn));
            gi2 = g2(uhi((Nn+1):2*Nn));
            gi = [gi1; gi2];

            ki = [k(uhi(1:Nn),uhi(Nn+1:2*Nn)); k(uhi(1:Nn),uhi(Nn+1:2*Nn))];
           
            vec_carg = A*uhi - dt*B*uhi + dt*C*gi-D*ki ; 

            uhi = A\vec_carg;

            T=0;
            N=0;



            % Representacion gráfica

                if mod(i,10) == 0

                    for j=1:length(t) 
                    
                    n1=t(1,j);
                    n2=t(2,j);
                    n3=t(3,j);
                
                    Atra=[p(1,n2)-p(1,n1) p(1,n3)-p(1,n1);p(2,n2)-p(2,n1) p(2,n3)-p(2,n1)];
                    F = @(x,y) Atra*[x;y]+[p(1,n1);p(1,n1)];
                    
                    w1=0.5*(1/3);
                    w2=w1;
                    w3=w1;
                    
                    p1 = F(0,0);
                    p2 = F(1,0);
                    p3 = F(0,1);
                    
                    T=T+det(Atra)*(w1*uhi(n1)+w2*uhi(n2)+w3*uhi(n3)); 
                    N=N+det(Atra)*(w1*uhi(Nn + n1)+w2*uhi(Nn + n2)+w3*uhi(Nn + n3)); 
                    
                    end                    

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

                     sgtitle(['tiempo = ', num2str(dt*i)]);

                     pause(0.001)
                end

end
