
%f(x) = 2x1^4 + 3x2^4 + 2x1^2 + 4x2^2 + x1*x2 -3*x1 -2*x2

%starting_point 
x = [0,0];

alfa = 0.1;
gamma = 0.9;
t_start = 1;

it = 1;
while not(isequal(gradient(x(1), x(2)),[0,0]))

    %risolvo il sistema lineare per trovare d
    d = -gradient(x(1),x(2))*inv(hessian(x(1),x(2)));

    t = t_start;
    
    vett = x + t*d;
    while f(vett(1), vett(2)) > f(x(1),x(2)) + alfa*t_start*d'*gradient(x(1), x(2))

        t = gamma*t;

    end

    x = x + t*d;

    fprintf("iteration %d with current point %s. Valore funzione: %s  valore gradiente: %s \n", it, mat2str(x), num2str(f(x(1),x(2))), gradient(x(1),x(2)));
    pause(1);
    
    it = it+1;


end



%ottieni il valore che l'Hessiana assume nel punto dato
function hessian_value = hessian(x1,x2)

    respect_x1 = 8*x1^3 + 4*x1 + x2 -3;
    respect_x2 = 12*x2^3 +8*x2 + x1 -2;

    respect_x1_x1 = 24*x1^2 + 4;
    respect_x1_x2 = 1;
    respect_x2_x1 = 1;
    respect_x2_x2 = 36*x2^2 + 8;

    hessian_value = [respect_x1_x1  respect_x1_x2; 
                     respect_x2_x1  respect_x2_x2];


end


%ritorna il valore della funzione
function f_value_in_x = f(x1,x2)

    f_value_in_x = 2*x1^4 + 3*x2^4 + 2*x1^2 + 4*x2^2 + x1*x2 -3*x1 -2*x2;

end


function gradient_in_x = gradient(x1,x2)

    respect_x1 = 8*x1^3 + 4*x1 + x2 -3;
    respect_x2 = 12*x2^3 +8*x2 + x1 -2;

    gradient_in_x = [respect_x1, respect_x2];

end