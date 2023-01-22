

%{
    min(x1 + 2x2 -3x3, -x1-x2-x3, -4x1-2x2+x3)
    x1 + x2 è x3 <= 10
    x3 <= 5
    x1,x2,x3 >=0

%}


%:::::::::::::::: TROVO IDEAL POINT ::::::::::::::::::::::::::::::::::::::

A = [1 1 1;
     0 0 1];

b =  [10; 5];

lb = [0 0 0];

%per ogni funzione trovo la sua soluzione ottima
z = [0;0;0];
for fi=1:3

    c = [1 2 -3];
    if fi==2
        c = [-1, -1 -1];
    elseif fi==3
        c = [-4 -3 1];
    end

    [xi,zi] = linprog(c,A,b,[],[],lb,[]);

    z(fi)= zi;
end


%:::::::::::::::::::::::: MOLP con norma 2 ::::::::::::::::::::::::::

C = [1   2 -3;
     -1 -1 -1;
     -4 -2  1];

H = C'*C;

c = -C'*z;

sol_2_norm = quadprog(H,c,A,b,[],[],lb,[]);


%:::::::::::::::::::::::: MOLP con norma 1 ::::::::::::::::::::::::::

c = [0 0 0 1 1 1]; %gli ultimi sarebbero gli y_i

A_2 = [C -eye(3)];
A_3 = [-C -eye(3)];
A_4 = [A_2; A_3; [A zeros(2,3)]];

b_1 = z;
b_2 = -z;
b_3 = [b_1; b_2; b];

sol_1_norm = linprog(c,A_4,b_3,[],[],lb,[]);


%:::::::::::::::::::::::: MOLP con norma infinito :::::::::::::::::::::::::

c = [0 0 0 1];

A_2 = [C -ones(3,1)];
A_3 = [-C -ones(3,1)];
A = [A_2; A_3; [A zeros(2,1)]];

sol_inf_norm = linprog(c,A,b_3,[],[],lb,[]);



%avendo ottenuto un risultato che sicuramente è weak per la inf norm voglio
%testare se tale risultato è anche di minimum risolvendo il problema
%ausiliario

f1 = @(x) x(1)+2*x(2)-3*x(3);
f2 = @(x) -x(1)-x(2)-x(3);
f3 = @(x) -4*x(1)-2*x(2)+x(3);

c = [0 0 0 1 1 1];

A = [1 1 1;
     0 0 1];

A_aux = [ [A zeros(2,3)]; 1 2 -3 1 0 0; -1 -1 -1 0 1 0; -4 -2 1 0 0 1];

b =  [10; 5];

b_aux = [b; f1(sol_inf_norm); f2(sol_inf_norm); f3(sol_inf_norm)];

lb = [ 0 0 0 0 0 0];

[sol, val] = linprog(c,A_aux,b_aux,[],[],lb,[]); %se 0 allora sol_inf_norm è di minimo e non solo weak

