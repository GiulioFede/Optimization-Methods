
%{
    Utilizzeremo una soft eps-SVM per trovare l'hyperpiano (quindi predittore
    lineare) che meglio predica i dati cercando di rispettare due
    condizioni tipiche dell'eps-SVM:
        1) massimo errore pari a eps=0.5
        2) predittore quanto più flat possibile per questione di
        generalizzazione.

    La differenza con la eps-SVM è che utilizziamo delle slack variable per
    permettere al sistema di trovare una soluzione anche quando eps è
    troppo piccolo.

    Qui i dati sono per non linearmente separabili e quindi utilizzeremo il
    kernel.
%}

dataset = [ 0    0.0620
    0.1000    0.2247
    0.2000    0.4340
    0.3000    0.6490
    0.4000    0.7370
    0.5000    0.8719
    0.6000    0.9804
    0.7000    1.0192
    0.8000    1.0794
    0.9000    1.0726
    1.0000    0.9252
    1.1000    0.8322
    1.2000    0.7457
    1.3000    0.5530
    1.4000    0.4324
    1.5000    0.2384
    1.6000    0.0060
    1.7000   -0.1695
    1.8000   -0.4023
    1.9000   -0.5487
    2.0000   -0.6583
    2.1000   -0.8156
    2.2000   -0.8582
    2.3000   -0.9217
    2.4000   -0.9478
    2.5000   -0.8950
    2.6000   -0.7947
    2.7000   -0.7529
    2.8000   -0.5917
    2.9000   -0.3654
    3.0000   -0.2392
    3.1000   -0.0172
    3.2000    0.2067
    3.3000    0.4111
    3.4000    0.5594
    3.5000    0.6678
    3.6000    0.7973
    3.7000    0.9605
    3.8000    1.0246
    3.9000    1.0947
    4.0000    1.0640
    4.1000    1.0070
    4.2000    0.9069
    4.3000    0.7604
    4.4000    0.6811
    4.5000    0.4661
    4.6000    0.2259
    4.7000    0.0944
    4.8000   -0.1224
    4.9000   -0.3606
    5.0000   -0.4550
    5.1000   -0.6669
    5.2000   -0.8049
    5.3000   -0.9114
    5.4000   -0.9498
    5.5000   -0.9771
    5.6000   -0.9140
    5.7000   -0.9127
    5.8000   -0.7953
    5.9000   -0.6653
    6.0000   -0.4486
    6.1000   -0.3138
    6.2000   -0.0900
    6.3000    0.0940
    6.4000    0.3098
    6.5000    0.4316
    6.6000    0.6899
    6.7000    0.8252
    6.8000    0.8642
    6.9000    0.9903
    7.0000    1.0232
    7.1000    1.0610
    7.2000    0.9887
    7.3000    0.9528
    7.4000    0.8486
    7.5000    0.7103
    7.6000    0.5312
    7.7000    0.3067
    7.8000    0.1591
    7.9000   -0.0511
    8.0000   -0.2771
    8.1000   -0.4264
    8.2000   -0.5930
    8.3000   -0.7232
    8.4000   -0.8070
    8.5000   -0.8913
    8.6000   -0.9097
    8.7000   -0.9874
    8.8000   -0.9269
    8.9000   -0.8212
    9.0000   -0.6551
    9.1000   -0.5258
    9.2000   -0.3894
    9.3000   -0.2136
    9.4000   -0.0436
    9.5000    0.2240
    9.6000    0.3940
    9.7000    0.5431
    9.8000    0.7247
    9.9000    0.8305
   10.0000    0.9881];

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
        M(row,column) = gaussian_kernel(x(row),x(column));   %<---------- NEW NEW NEW
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




%trovo quel lambda (positivo) che sia compreso tra 0 e C (entrambi esclusi)
index = find(lambda(1:l)>1e-4 & lambda(1:l)<C-1e-4);
b=0;
%se non è vuoto...
if length(index)>0
    index = index(1);
    %non ho w quindi calcolo k(x_j, x_i) dove x_i è fisso
    k_values = arrayfun(@(index_j) gaussian_kernel(x(index), x(index_j)), 1:l);
    b = y(index)-(sum( (lambda(1:l)-lambda(l+1:end)).*k_values'))-eps;
%...altrimenti utilizzo i lambda negativi
else
    index = find(lambda(l+1:end)>1e-4 & lambda(l+1:end)<C-1e-4);
    index = index(1);
    %non ho w quindi calcolo k(x_j, x_i) dove x_i è fisso
    k_values = arrayfun(@(index_j) gaussian_kernel(x(index), x(index_j)), 1:l);
    b = y(index)-(sum( (lambda(1:l)-lambda(l+1:end)).*k_values'))-eps;
end


%plotto. 
start_p = min(x)-5;
end_p = max(x)+5;
points = start_p:0.1:end_p;
y_points = zeros(length(points),1);

for p=1:length(points)
    
    %calcolo la predizione sul punto
    k_values = arrayfun(@(index_j) gaussian_kernel(x(index_j), points(p)), 1:l);
    %predico y
    y_points(p) = (sum( (lambda(1:l)-lambda(l+1:end)).*k_values')) + b;
end


hold on
plot([points], [y_points],'color','r');

%evidenzio i support vector. Sono quelli con uno dei due lambda relativi positivi
index_lambda_pos = find(lambda(1:l)>1e-3);
index_lambda_neg = find(lambda(l+1:end)>1e-3);
%prendo l'unione
common = union(index_lambda_pos, index_lambda_neg);

scatter(x(common), y(common), 'green');



function k = gaussian_kernel(x1, x2)
    gamma = 1;
    k = exp(-gamma*norm(x1-x2)^2);

end

function k = polynomial_kernel(x1, x2)
    p = 4;
    k = (x1'*x2+1)^p;

end
