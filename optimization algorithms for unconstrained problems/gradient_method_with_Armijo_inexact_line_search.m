

%f(x) = 2x1^4 + 3x2^4 + 2x1^2 + 4x2^2 + x1*x2 -3*x1 -2*x2

% starting point
x = [0 0];

alfa = 0.1;
gamma = 0.9;

fprintf("Starting point %s. Valore funzione: %s \n", mat2str(x), num2str(f(x(1), x(2))));

it = 1;
old_f_value = inf;
while true

    %calcolo vettore gradiente in x corrente
    gradient_in_x = gradient(x(1), x(2));

    %prendo come direzione quella opposta
    direction = -gradient_in_x;

    %calcolo un approssimazione dello step size con la Armijo inexact line
    %search
    t = 1;
    vett1 = x + t*direction;
    while (f(vett1(1), vett1(2)) > ( f(x(1),x(2)) + alfa*t*direction'*gradient(x(1),x(2)) )  )
        
        t = gamma*t;

        vett1 = x + t*direction;
        
    end
    step = t;
    
    pause(0.3)

    %calcolo il nuovo punto
    x = x + step*direction;

    current_f_value = f(x(1), x(2));
    fprintf("iteration %d with current point %s. Valore funzione: %s \n", it, mat2str(x), num2str(current_f_value));
    it = it+1;


    if abs(current_f_value-old_f_value)<1e-6
        break;
    else
        old_f_value = current_f_value;
    end
end


%ritorna il valore della funzione
function f_value_in_x = f(x1,x2)

f_value_in_x = 2*x1^4 + 3*x2^4 + 2*x1^2 + 4*x2^2 + x1*x2 -3*x1 -2*x2;

end

%Calcolo il vettore gradiente come una funzione
function gradient_in_x = gradient(x1,x2)

    respect_x1 = 8*x1^3 + 4*x1 + x2 -3;
    respect_x2 = 12*x2^3 +8*x2 + x1 -2;

    gradient_in_x = [respect_x1, respect_x2];

end