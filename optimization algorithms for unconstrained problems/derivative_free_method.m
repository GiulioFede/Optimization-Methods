
%f(x) = 2x1^4 + 3x2^4 + 2x1^2 + 4x2^2 + x1*x2 -3*x1 -2*x2

%starting point
x = [0,0];
t = 5;
beta = 0.5;
epsilon = 1e-5;
base = [ 1 0; 0 1; -1 -1];
number_of_positive_vectors = size(base,1);

it = 1;
while t>epsilon

    %genero i punti su ogni vettore di base
    pool_set = zeros(number_of_positive_vectors,2);
    %valuto ciascun punto del pool set
    pool_set_values = zeros(number_of_positive_vectors,1);

    for i=1:number_of_positive_vectors
        pool_point = x + t*base(i,:);
        pool_set(i,:) = pool_point;
        pool_set_values(i) = f(pool_point(1), pool_point(2));
    end

    %cerco il punto che ha valore minore di quello attuale
    found = false;
    for i=1:number_of_positive_vectors
        if pool_set_values(i) < f(x(1), x(2))
            found = true;
            x = x + t*base(i,:);
            break;
        end
    end

    %se non ho trovato nessuno
    if found==false
        t = beta*t;
    end

    fprintf("iteration %d with current point %s. Valore funzione: %s  valore step: %s \n", it, mat2str(x), num2str(f(x(1),x(2))), num2str(t));
    
    it = it+1;

    pause(1);

end



%ritorna il valore della funzione
function f_value_in_x = f(x1,x2)

    f_value_in_x = 2*x1^4 + 3*x2^4 + 2*x1^2 + 4*x2^2 + x1*x2 -3*x1 -2*x2;

end