% --- Interpolación funciones R1 -> R1 con elementos finitos de orden m uniformemente distribuidos---


f = @(x) 1./(1+25*x.^2); %función a interpolar

a=-1; 
b=1;

Ne = 10; %numero de elementos
m=1; %orden de la interpolación polinómica

h = (b-a)/Ne;
ai = a;
bi= a+h;

for k= 1:Ne

    x = ai:0.001:bi; %puntos en los que obtener evaluado el pol interpolador (en cada elemento)
    p = 0*x; 

    xi=ai:(bi-ai)/m:bi;
    fi = f(xi);

    for i=1:length(xi)
        li=0*x+1;

        for j=1:length(xi)
            if i ~= j
                li = li.*(x-xi(j))/(xi(i)-xi(j));
            end
        end

        p = p + fi(i).*li;

    end

    plot(x,p, 'b') %plot interpolación
    %scatter(xi,fi)
    hold on

    ai = ai + h;
    bi = bi + h;

end

z = a:0.001:b;
plot(z,f(z),'g') %plot función 

