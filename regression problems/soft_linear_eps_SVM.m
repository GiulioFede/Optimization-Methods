
%{
    Utilizzeremo una soft eps-SVM per trovare l'hyperpiano (quindi predittore
    lineare) che meglio predica i dati cercando di rispettare due
    condizioni tipiche dell'eps-SVM:
        1) massimo errore pari a eps=0.5
        2) predittore quanto più flat possibile per questione di
        generalizzazione.

    La differenza con la eps-SVM è che utilizziamo delle slack variable per
    permettere al sistema di trovare una soluzione anche quando eps è
    troppo piccolo
%}

dataset = [
     0    2.5584
    0.5000    2.6882
    1.0000    2.9627
    1.5000    3.2608
    2.0000    3.6235
    2.5000    3.9376
    3.0000    4.0383
    3.5000    4.1570
    4.0000    4.8498
    4.5000    4.6561
    5.0000    4.5119
    5.5000    4.8346
    6.0000    5.6039
    6.5000    5.5890
    7.0000    6.1914
    7.5000    5.8966
    8.0000    6.3866
    8.5000    6.6909
    9.0000    6.5224
    9.5000    7.1803
   10.0000    7.2537];

%plotto i punti
scatter(dataset(:,1), dataset(:,2), 'black');

%dimensione del dataset
l = size(dataset,1);

%osservazioni
x = dataset(:,1);

%ground truth
y= dataset(:,2);

%fisso l'epsilon e lo metto davvero piccolo
eps = 0.2;

C = 10;

% la quadprog vuole come funzione obiettivo 1/2 x^T H x + f^Tx
H = [1 0 zeros(1,2*l); 0 0 zeros(1,2*l); zeros(2*l,2*l+2)];
f = [zeros(2,1); C*ones(2*l,1)];

I = eye(l);
%i vincoli desiderati sono A x <= b
A = [-x -ones(l,1) -I zeros(l); x ones(l,1) zeros(l) -I ];
b = [-y+eps; y+eps];

%lower bound
lb = [-inf, -inf, zeros(1,2*l)];

%risolvo. La soluzione è il vettore z=(w,b)
z = quadprog(H,f,A,b,[],[],lb,[]);
w = z(1);
b = z(2);

%plotto. 
start_p = min(x)-5;
end_p = max(x)+5;
y_start = w*start_p + b;
y_end = w*end_p + b;

hold on
plot([start_p, end_p], [y_start, y_end],'color','r');
