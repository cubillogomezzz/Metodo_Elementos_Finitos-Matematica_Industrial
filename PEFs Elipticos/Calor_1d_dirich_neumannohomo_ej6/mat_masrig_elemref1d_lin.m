% --- Cálculo matrices masas y rigideces en elemento de referencia 1D mediante elem lineales --- 

% Utilizando cuadratura gauss-leg

phi0 = @(x) 0.5*(1-x);
phi1 = @(x) 0.5*(1+x);

dphi0 = @(x) 0*x-0.5; 
dphi1 = @(x) 0*x + 0.5;

xgi = [-sqrt(3)/3 sqrt(3)/3]; % Nodos y pesos GL
wi = [1 1];

Mg = zeros(2,2);
Rg = Mg;

phi0i = phi0(xgi);
phi1i = phi1(xgi);

dphi0i = dphi0(xgi);
dphi1i = dphi1(xgi);

Mg(1,1) = sum(wi .* phi0i .* phi0i);
Mg(1,2) = sum(wi .* phi0i .* phi1i); 
Mg(2,2) = sum(wi .* phi1i .* phi1i);
Mg(2,1) = Mg(1,2);

Rg(1,1) = sum(wi .* dphi0i .* dphi0i);
Rg(1,2) = sum(wi .* dphi0i .* dphi1i); 
Rg(2,2) = sum(wi .* dphi1i .* dphi1i);
Rg(2,1) = Rg(1,2);
