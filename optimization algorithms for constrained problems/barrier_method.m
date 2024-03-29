%{
        Problema:
                min  1/2 (x1-3)^2 + (x2-2)^2
                     
                      -2x1 + x2 <=0
                      x1 + x2 <= 4
                      -x2 <= 0

%}

%{
    ho bisogno di un punto interno alla feasible region, quindi risolvo il
    problema ausiliario:
                    min s
                    g_i(x)<=s
%}
tolerance = 1e-6;
tau = 0.5;
global epsilon;

epsilon = 1;
%numero vincoli
m = 3;

%prendo un qualsiasi punto iniziale, come (1,1)
x = [6,7];
%lo valuto su tutti i vincoli
vv1 = -2*x(1) + x(2);
vv2 = x(1) + x(2) - 4;
vv3 = -x(2);
%prendo il max
maximum = max([vv1,vv2,vv3]);
%prendo s maggiore di tale maximum
s = maximum + 0.5;

%utilizzo (x,s) come punto interno alla feasible region del problema
%ausiliario
x = [x s];

%utilizo barrier method
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
k=1;
while true

    %trovo l'ottima soluzione del problema
    f_and_g = @get_function_and_gradient_of_auxiliary_problem;
    
    x = fminunc(f_and_g, x, options);
    fprintf("soluzione iterazione %d: %s   ", k, mat2str(x));
    
    if m*epsilon < tolerance
        fprintf("\n \nTerminato. Soluzione finale: %s   \n\n", mat2str(x));
        break;
    else
        k = k+1;
        epsilon = tau*epsilon;
    end
end

%controlliamo se s<0, ossia esiste l'interior
if x(3)<0
    fprintf("Interior esistente: %s   \n", mat2str(x(1:2)));

    %quindi utilizzo le prime due componenti come starting point
    x = x(1:2);

    %resetto
    epsilon = 1;
    k = 1;

    while true
        %trovo l'ottima soluzione del problema
        f_and_g = @get_function_and_gradient;
        
        x = fminunc(f_and_g, x, options);
        fprintf("soluzione iterazione %d: %s   ", k, mat2str(x));
        
        if m*epsilon < tolerance
            fprintf("\n \nTerminato. Soluzione finale: %s   \n\n", mat2str(x));
            break;
        else
            k = k+1;
            epsilon = tau*epsilon;
        end

    end


else
    fprintf("Interior non esistente. \n");
end



function [objective_function,gradient] = get_function_and_gradient_of_auxiliary_problem(x)

    global epsilon;

    % f = s - epsilon*sum(log(-g_i)  ricorda però che qui g_i(x)<=s e non <=0
    objective_function = x(3) - epsilon*( log(-(-2*x(1) + x(2)-x(3))) + log(-(x(1) + x(2) - 4-x(3))) +  log(-(-x(2)-x(3))));
    
    %gradient
    gradient = [0,0,1] + epsilon*( (1/-(-2*x(1) + x(2)-x(3)))*[-2,1,-1] + (1/-(x(1) + x(2) - 4-x(3)))*[1,1,-1] + (1/-(-x(2)-x(3)))*[0,-1,-1] );

    
end


function [objective_function,gradient] = get_function_and_gradient(x)

    global epsilon;

    % f = s - epsilon*sum(log(-g_i)  ricorda però che qui g_i(x)<=s e non <=0
    objective_function = 1/2*(x(1)-3)^2 + (x(2)-2)^2 - epsilon*( log(-(-2*x(1) + x(2))) + log(-(x(1) + x(2) - 4)) +  log(-(-x(2))));
    
    %gradient
    gradient = [x(1)-3,2*x(2)-4] + epsilon*( (1/-(-2*x(1) + x(2)))*[-2,1] + (1/-(x(1) + x(2) - 4))*[1,1] + (1/-(-x(2)))*[0,-1] );

    
end
