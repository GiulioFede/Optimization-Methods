
%{
        Problema:
                min  1/2 (x1-3)^2 + (x2-2)^2
                     
                      -2x1 + x2 <=0
                      x1 + x2 <= 4
                      -x2 <= 0

%}

global epsilon;
epsilon = 5;
tao = 0.5;

k=1;


while true

    %{

        Adesso devo risolvere il problema unconstrained. Utilizzerò una
        funzione di Matlab che fondamentalmente fa quanto fatto nella
        cartella "optimization algorithms for unconstrained problems",
        quindi volendo potremmo anche utilizzare gli algoritmi scritti
        prima.

        Tale funzione vuole:
            1) funzione obiettivo da minimizzare su tutto lo spazio
            2) (facoltativo MA importante) il gradiente della funzione
            obiettivo. Anche se facoltativo lo forniamo dato che molti
            metodi sono basati sul gradiente quindi la funzione proverà a
            creare una approssimazione del gradiente della funzione
            obiettivo, ovviamente avrà meno precisione rispetto a una forma
            esplicita. Per questo la daremo.
        
        NB: devo creare una funzione matlab che restituisca f e g con f
        funzione obiettivo e g gradiente relativo (la trovi alla fine)
    %}

    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
    
    f_and_g = @get_function_and_gradient_of_penalized_problem;
    %risolvo il penalty problem
    x = fminunc(f_and_g, x, options);
    fprintf("soluzione iterazione %d: %s   ", k, mat2str(x));

    %controllo sel la soluzione rispetta i vincoli, quindi appartiene alla
    %feasible region, perchè se è cosi allora tale soluzione è anche
    %soluzione del problema originale
    
    %per ogni vincolo
    [check, max_residual] = check_if_all_constraints_are_respected(x(1), x(2));

    if(check == true)
        fprintf("Soluzione ottima trovata: %s\n.", mat2str(x));
        break;
    else
        fprintf("La soluzione trovata non appartiene alla feasible region: epsilon %d,  max residual: %d.\n", epsilon, max_residual);
    end

    k= k+1;
    epsilon = tao*epsilon;

    %se comunque max(Ax_k -b) è piccolo, mi fermo
    if max_residual < 1e-3
        fprintf("Differenza tra Ax e b sotto la soglia. Soluzione finale: %s.\n", mat2str(x));
        break;
    end

end



function [objective_function,gradient] = get_function_and_gradient_of_penalized_problem(x)

    global epsilon;

    %trovo la soluzione del penalized problem  f(x) + 1/epsilon p(x) dove
    %p(x) = sum( max(0,g_i(x)^2)

    penalty_function = 1/epsilon * (max(0, -2*x(1)+x(2))^2 + max(0,x(1)+x(2)-4)^2 + max(0,-x(2))^2);

    %funzione obiettivo del nuovo penalized problem
    objective_function = 1/2*(x(1)-3)^2 + (x(2)-2)^2  +  penalty_function;

    %trovo gradiente (dato che f (originale) e i vincoli g_i sono
    %continuamente differenziabili allora il gradiente della funzione
    %obiettivo con penalty function ha una forma particolare
    gradient_of_f_original = [ x(1)-3   2*x(2)-4];

    gradient = gradient_of_f_original + 1/epsilon * (max(0, -2*x(1)+x(2))*[-2 1] + max(0,x(1)+x(2)-4)*[1 1] + max(0,-x(2))*[0 -1]);

end

function [check,max_residual] = check_if_all_constraints_are_respected(x1,x2)

    A = [-2 1; 1 1; 0 -1];
    b = [0 4 0];

    b1 = A*[x1 x2]';
    max_residual = max(b1'-b);
    if (sum(b1<=b)==3)
        check = true;
    else
        check = false;
    end

end