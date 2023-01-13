
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

%{
    la quadprog vuole come funzione obiettivo 1/2 l^T H l + c^T l
    Ricordo che: 
                Q= [ M    -M
                    -M     M ]

    mentre:
                c = [(-eps+y_1),...,(eps-y_l)]
%}
M = zeros(l,l);
for row=1:l
    for column=1:l
        M(row,column) = x(row)*x(column);
    end
end
Q = [M  -M;  -M  M];

%c = [-1*ones(l,1)*eps+y; ones(l,1)*eps-y];
c = [-1*ones(l,1)*eps+y; -ones(l,1)*eps-y];


%abbiamo un solo vincolo di uguaglianza --> Aeq*l=beq
Aeq= [ones(1,l)  -1*ones(1,l)];
beq = 0;

%lower bound
lb = zeros(2*l,1);
up = C*ones(2*l,1);

%risolvo. La soluzione sono i vettori lambda (2*l) positivi e negativi
%ricorda di cambiare il segno della funzione obiettivo in quanto il
%problema duale è espresso come "max".
lambda = quadprog(Q,-c,[],[],Aeq,beq,lb,up);


w = sum( (lambda(1:l)-lambda(l+1:end)).*x);

%trovo quel lambda (positivo) che sia compreso tra 0 e C (entrambi esclusi)
index = find(lambda(1:l)>=0.001 & lambda(1:l)<C-1e-3);
%se non è vuoto...
if length(index)>0
    index = index(1);
    b = y(index)-w*x(index)-eps;
%...altrimenti utilizzo i lambda negativi
else
    index = find(lambda(l+1:end)>=0.001 & lambda(l+1:end)<C-1e-3);
    index = index(1);
    b = y(index)-w*x(index)-eps;
end


%plotto. 
start_p = min(x)-5;
end_p = max(x)+5;
y_start = w*start_p + b;
y_end = w*end_p + b;

hold on
plot([start_p, end_p], [y_start, y_end],'color','r');

%evidenzio i support vector. Sono quelli con uno dei due lambda relativi positivi
index_lambda_pos = find(lambda(1:l)>1e-3);
index_lambda_neg = find(lambda(l+1:end)>1e-3);
%prendo l'unione
common = union(index_lambda_pos, index_lambda_neg);

scatter(x(common), y(common), 'green');
