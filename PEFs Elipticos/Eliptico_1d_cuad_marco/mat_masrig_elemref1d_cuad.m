% --- Cálculo matrices masas y rigideces en elemento de referencia 1D mediante elem cuadráticos --- 

% Utilizando cuadratura gauss-leg

phi0 = @(x) 0.5*(x-1).*x;
phi1 = @(x) 1-x.^2;
phi2 = @(x) 0.5*(1+x).*x;

dphi0 = @(x) x-0.5; 
dphi1 = @(x) -2*x; 
dphi2 = @(x) x + 0.5;

xgi = [-sqrt(3)/3 sqrt(3)/3]; % Nodos y pesos GL
wi = [1 1];

Mg = zeros(3,3);
Rg = Mg;

phi0i = phi0(xgi);
phi1i = phi1(xgi);
phi2i = phi2(xgi);

dphi0i = dphi0(xgi);
dphi1i = dphi1(xgi);
dphi2i = dphi2(xgi);

Mg(1,1) = sum(wi .* phi0i .* phi0i);
Mg(2,2) = sum(wi .* phi1i .* phi1i);
Mg(3,3) = sum(wi .* phi2i .* phi2i);

Mg(1,2) = sum(wi .* phi0i .* phi1i); 
Mg(1,3) = sum(wi .* phi0i .* phi2i);
Mg(2,3) = sum(wi .* phi1i .* phi2i);

Mg(2,1) = Mg(1,2);
Mg(3,1) = Mg(1,3);
Mg(3,2) = Mg(2,3);

Rg(1,1) = sum(wi .* dphi0i .* dphi0i);
Rg(2,2) = sum(wi .* dphi1i .* dphi1i);
Rg(3,3) = sum(wi .* dphi2i .* dphi2i);

Rg(1,2) = sum(wi .* dphi0i .* dphi1i);
Rg(1,3) = sum(wi .* dphi0i .* dphi2i);
Rg(2,3) = sum(wi .* dphi1i .* dphi2i);

Rg(2,1) = Rg(1,2);
Rg(3,1) = Rg(1,3);
Rg(3,2) = Rg(2,3);



