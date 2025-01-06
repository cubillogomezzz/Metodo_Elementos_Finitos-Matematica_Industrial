% --- Cálculo matrices de rigidez, masas y convección ---

% La función [R M b] = assema(p,t,a,k,f) calcula las matrices de
% rigiez, masas, y el vector de carga para un problema eliptico puramente
% difusivo (ec. del calor) con k matricial. En nuestra implementación no 
% utilizamos assema para calcular b, el cual obtenemos externamente a 
% partir de R y M. Tampoco incluimos a y k, las cuales fijamos en 1 y 
% tratamos externamente. 

% Matlab no tiene una forma cómoda de calcular matrices de convección. Nos
% aprovechamos de assema utilizando el output del vector de carga, que se
% calcula como bj = integral(f*phi_j). Si en f introducimos  
% parcial_x1(phi_i), obtenemos Cx_ji. En lugar de usar un function handle,
% sabiendo que assema integra con una gauss-legendre de un punto en el 
% centro de cada triángulo, sirve introducir el valor de parcial_x1(phi_i) 
% en dicho punto de integración. La función pdecgrad(p,t,u,c) hace esto.

Cx = sparse(length(xi),length(xi));
Cy = Cx;

for i = 1:length(xi)
    phi_i = 0*xi';
    phi_i(i) = 1;
    [phi_i_x phi_i_y] = pdecgrad(p,t,'1',phi_i); %¿como funciona pdecgrad?
    [R M b_gx] = assema(p,t,1,1,phi_i_x);
    Cx(:,i) = b_gx;
    [R M b_gy] = assema(p,t,1,1,phi_i_y);
    Cy(:,i) = b_gy;
end
