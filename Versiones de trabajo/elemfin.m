% --- Interpolación Elem Fin ---

%equiespaciado todo

f = @(x) 1./(1+25*x.^2); %función conocida o a comparar con la interpolación
% fi = [ ];


a=-1;
b=1;




Ne = 10; %param entrada
h = (b-a)/Ne;
m=3; 

ai = a;
bi= a+h;

for k= 1:Ne

    xi=ai:(bi-ai)/m:bi;
    fi = f(xi);
    x = ai:0.001:bi; %puntos en los que obtener evaluado el pol interpolador
    p = 0*x; 

                for i=1:length(xi) 
                    li=0*x+1;
                
                    for j=1:length(xi)
                        if i ~= j
                         li = li.*(x-xi(j))/(xi(i)-xi(j));
                        end
                    end
                
                    p = p + fi(i).*li; 
                
                end

                plot(x,p, 'b')
                plot(xi,fi)
                hold on

                ai = ai +h;
                bi = bi +h;

end
 
