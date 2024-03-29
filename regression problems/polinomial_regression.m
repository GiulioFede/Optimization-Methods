%{

    Risolveremo il problema di regressione utilizzando un polinomio di 
    grado 3.
    Inoltre lo risolveremo nei 3 casi in cui si utilizzi la norma 1, 2 e la
    infinity.

%}

dataset = [
    -5.0000  -96.2607 
   -4.8000  -85.9893
   -4.6000  -55.2451
   -4.4000  -55.6153
   -4.2000  -44.8827
   -4.0000  -24.1306
   -3.8000  -19.4970
   -3.6000  -10.3972
   -3.4000   -2.2633
   -3.2000    0.2196
   -3.0000    4.5852
   -2.8000    7.1974
   -2.6000    8.2207
   -2.4000   16.0614
   -2.2000   16.4224
   -2.0000   17.5381
   -1.8000   11.4895
   -1.6000   14.1065
   -1.4000   16.8220
   -1.2000   16.1584
   -1.0000   11.6846
   -0.8000    5.9991
   -0.6000    7.8277
   -0.4000    2.8236
   -0.2000    2.7129
         0    1.1669
    0.2000   -1.4223
    0.4000   -3.8489
    0.6000   -4.7101
    0.8000   -8.1538
    1.0000   -7.3364
    1.2000  -13.6464
    1.4000  -15.2607
    1.6000  -14.8747
    1.8000   -9.9271
    2.0000  -10.5022
    2.2000   -7.7297
    2.4000  -11.7859
    2.6000  -10.2662
    2.8000   -7.1364
    3.0000   -2.1166
    3.2000    1.9428
    3.4000    4.0905
    3.6000   16.3151
    3.8000   16.9854
    4.0000   17.6418
    4.2000   46.3117
    4.4000   53.2609
    4.6000   72.3538
    4.8000   49.9166
    5.0000   89.1652];

%plotto i punti
scatter(dataset(:,1), dataset(:,2), 'black');

%dimensione del dataset
l = size(dataset,1);

%osservazioni
x = dataset(:,1);

%ground truth
y= dataset(:,2);

%fisso il grado del polinomio
degree_of_polynomial = 3;


%---------------------------- Norma 2 ------------------------------------
%{
    Nel caso della norma 2, il problema (least squares problem) è risolto
    risolvendo il sistema lineare:
                            A^T A z = A^T y

    La matrice A ha per ogni riga i il punto i e in ogni
    colonna j tale punto è elevato alla j-1. Mentre il vettore y è
    quello dei valori che ogni punto deve assumere.
%}

A = [ones(l,1),x,x.^2, x.^3]; %perchè voglio un polinomio di grado 3



%risolvo il sistema A^T A z = A^T y con linsolve
z = linsolve(A'*A, A'*y);

%plotto il polinomio che fitta i dati
min_x = min(x)-5;
max_x = max(x)+5;
points = min_x:0.01:max_x;
predictions = arrayfun(@(p) (z(1) + z(2)*p + z(3)*p^2 + z(4)*p^3), points);

hold on
plot(points, predictions,'color','r');


%---------------------------- Norma 1 ------------------------------------
%{
    Nel caso della norma 1, il problema diventa lineare (dopo vari trick 
    applicati al sistema di partenza).
                            min u1+u2+...+ul  + 0*z1+ 0*z2 + ... +0*zl
                            ui >= Ai z - yi
                            ui >= yi - Ai z

    Utilizzeremo linprog il quale vuole:
    1) f:
        la funzione obiettivo è espressa come f^T x. Nel nostro caso f è un
        vettore tutti di 1. Nota che f include anche gli zeri per gli z
    2) A x <= b:
        Compatto tutto come M <= b in questo modo:
                        [ A  -I ] [z]       [y]
                                        <= 
                        [ -A -I ] [u]       [-y]
%}

f = [zeros(4,1);ones(l,1)];

I = eye(l);
M = [ [A; -A], [-I; -I]];

b = [y;-y];

sol = linprog(f,M,b);

%prendo solo le prime 4 componenti, quelle relative a z
z = sol(1:4);

%plotto
predictions = arrayfun(@(p) (z(1) + z(2)*p + z(3)*p^2 + z(4)*p^3), points);

hold on
plot(points, predictions,'color','b');



%---------------------------- Norma infinity ------------------------------------
%{
    Stessa cosa della norma 1, solo che qui abbiamo solo una u, vedila come
    una u con una sola componente
%}

f = [zeros(4,1); 1];

M = [ [A; -A], -ones(2*l,1)];

b = [y;-y];

sol = linprog(f,M,b);

%prendo solo le prime 4 componenti, quelle relative a z
z = sol(1:4);

%plotto
predictions = arrayfun(@(p) (z(1) + z(2)*p + z(3)*p^2 + z(4)*p^3), points);

hold on
plot(points, predictions,'color','g');



legend('observations','norm 1','norm 2', 'norm infinity');

