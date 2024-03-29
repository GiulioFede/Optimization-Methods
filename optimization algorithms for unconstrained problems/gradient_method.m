

%f(x) = 3x1^2 +3x2^2 +3x3^2 +3x4^2-4x1x3-4x2x4 + x1 - x2 +2x3 -3x4

%se è a true utilizziamo il conjugate gradient method
use_conjugate = true;

% starting point
x = [0,0,0,0];

%scrivo l'hessiana Q
Q = [6 0 -4 0;
     0 6  0 -4;
     -4 0 6 0;
     0 -4 0  6];

fprintf("Starting point %s. Valore funzione: %s \n", mat2str(x), num2str(f(x(1), x(2), x(3), x(4))));

it = 1;
old_f_value = inf;
prev_direction = -1; %utile solo per il conjugate 
while true

    %calcolo vettore gradiente in x corrente
    gradient_in_x = gradient(x(1), x(2), x(3), x(4));

    %se uso il conjugate method
    if use_conjugate == true && it>1
        beta = (gradient_in_x*Q*prev_direction')/(prev_direction*Q*prev_direction');
        direction = -gradient_in_x + beta*prev_direction;
    else
        %altrimenti...
        %prendo come direzione quella opposta
        direction = -gradient_in_x;
    end

    %calcolo il migliore step size
    step = -(gradient_in_x * direction')/(direction*Q*direction');

    %calcolo il nuovo punto
    x = x + step*direction;

    current_f_value = f(x(1), x(2), x(3), x(4));
    fprintf("iteration %d with current point %s. Valore funzione: %s \n", it, mat2str(x), num2str(current_f_value));
    it = it+1;

    prev_direction = direction;


    if abs(current_f_value-old_f_value)<1e-6
        break;
    else
        old_f_value = current_f_value;
    end

    pause(0.3);
end


%ritorna il valore della funzione
function f_value_in_x = f(x1,x2,x3,x4)

f_value_in_x = 3*x1^2 +3*x2^2 +3*x3^2 +3*x4^2-4*x1*x3-4*x2*x4 + x1 - x2 +2*x3 -3*x4;

end

%Calcolo il vettore gradiente come una funzione
function gradient_in_x = gradient(x1,x2,x3,x4)

    respect_x1 = 6*x1 -4*x3 +1;
    respect_x2 = 6*x2 -4*x4 -1;
    respect_x3 = 6*x3 -4*x1 +2;
    respect_x4 = 6*x4 -4*x2 -3;

    gradient_in_x = [respect_x1, respect_x2, respect_x3, respect_x4];

end