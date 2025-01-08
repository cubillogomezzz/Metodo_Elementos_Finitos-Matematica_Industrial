%------ PROBLEMA PARABÓLICO -----------%
%
% Tt = Dt*laplaciano(T) + at*(1-T/kt)*T - alphatn*N*T
% Nt = Dn*laplaciano(N) + an*(1-N/kn)*N - alphant*N*T
% dT/dn(fron) = 0
% dN/dn(fron) = 0
% T(x,0) = T0
% N(x,0) = N0
%
%------ discretización: crank-nicolson ---------%

%load('mallado_tumor.mat')
an = 5;   %[1/día]
at = 0;    %[1/día]
Dn = 2.0e-15;   %[cm^2/día]
Dt = 4.2e-3;    %[cm^2/día]
kn = 1; %capacidad carga del medio
kt = 1;
alpha_tn = 1.0e-4;   %[1/día]
alpha_nt = 0.2;   %[1/día]

%parametros en el tiempo
dt = 0.01;
t0 = 0;
tf = 10;
Nt = (tf - t0) / dt;

%parametros de newton-raphson
N_newton = 10000;
delta = 1e-10;

%para el mallado usamos pdetool
%load('mallado_tumor.mat')
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)'; %transponemos para usar el comando trisurf
%frontera de D
fron = unique([e(1,:) e(2,:)]);
Nn = length(xi); %numero de nodos

%assema calcula las matrices de rigidez y de masas
[R M] = assema(p,t,1,1,0);

%jacobiana  del sistema para newton-raphson
I = eye(Nn);
JF = @(X) [(1-0.5*dt*at)*M+0.5*dt*Dt*R,0*M,0.5*dt*at/kt*M,0*M,0.5*dt*alpha_tn*M;...
    0*M,(1-0.5*dt*an)*M+0.5*dt*Dn*R,0*M,0.5*dt*an/kn*M,0.5*dt*alpha_nt*M;...
    -2*diag(X(1:Nn)),0*M,I,0*M,0*M;...
    0*M,-2*diag(X(Nn+1:2*Nn)),0*M,I,0*M;...
    -diag(X(Nn+1:2*Nn)),-diag(X(1:Nn)),0*M,0*M,I];

%condiciones iniciales de las variables
T0 = @(x,y) exp(-((x-0.5).^2+(y-0.5).^2)/0.01);
Thn = T0(xi,yi)';
Nhn = 0*xi'+1;
whn = Thn .* Nhn;
z_thn = Thn.^2;
z_nhn = Nhn.^2;
figure(1)
trisurf(elem,xi,yi,Thn);
title(0)
% view(2)
% shading interp
% colorbar
% caxis([-0.1 0.1])
axis([-1 1 -1 1 0 1])
pause(1e-10)
for n = 1:Nt

    F = @(X) [M*((1-0.5*dt*at)*X(1:Nn)+0.5*dt*at/kt*X(2*Nn+1:3*Nn)+0.5*dt*alpha_tn*X(4*Nn+1:5*Nn))+R*X(1:Nn)*0.5*dt*Dt-M*((1+0.5*dt*at)*Thn-0.5*dt*at/kt*z_thn-0.5*dt*alpha_tn*whn)-R*Thn*0.5*dt*Dt;...
        M*((1-0.5*dt*an)*X(Nn+1:2*Nn)+0.5*dt*an/kn*X(3*Nn+1:4*Nn)+0.5*dt*alpha_nt*X(4*Nn+1:5*Nn))+R*X(Nn+1:2*Nn)*0.5*dt*Dn-M*((1+0.5*dt*an)*Nhn-0.5*dt*an/kn*z_nhn-0.5*dt*alpha_nt*whn)-R*Nhn*0.5*dt*Dn;...
        X(2*Nn+1:3*Nn)-X(1:Nn).^2;...
        X(3*Nn+1:4*Nn)-X(Nn+1:2*Nn).^2;...
        X(4*Nn+1:5*Nn)-X(1:Nn).*X(Nn+1:2*Nn)];
    Xn = [Thn;Nhn;z_thn;z_nhn;whn];
    for i = 1:N_newton
        Xn_ant = Xn;
        Xn = Xn - JF(Xn)\F(Xn);
        if norm(Xn-Xn_ant)/norm(Xn)<delta
            break;
        end
    end
    Thn = Xn(1:Nn);
    Nhn = Xn(Nn+1:2*Nn);
    z_thn = Xn(2*Nn+1:3*Nn);
    z_nhn = Xn(3*Nn+1:4*Nn);
    whn = Xn(4*Nn+1:5*Nn);

    if mod(n,10) == 0
        figure(1)
        subplot(2, 1, 1)
        trisurf(elem,xi,yi,Thn)
        title('CRANK NICHOLSON',n*dt)
        shading interp
        colorbar
        caxis([0 1])
        axis([0 1 0 1 0 1])
        

        subplot(2, 1, 2)
        trisurf(elem,xi,yi,Nhn)
        title('CRANK NICHOLSON',n*dt)
        shading interp
        colorbar
        caxis([0 1])
        axis([0 1 0 1 0 1])
        pause(0.01)
    end

    
end
