
%{
    Questo algoritmo prevede lo stesso obiettivo fatto nell'omonimo,
    ma stavolta risolvendo il più semplice problema Lagrangiano duale.

    Inoltre il dataset contiene dati NON linearmente separabili.
    Utilizzando il duale Lagrangiano c'è poco da modificare. 

%}




%{
    DATASET
%}

A = [0.4952    6.8088
        2.6505    8.9590
        3.4403    5.3366
        3.4010    6.1624
        5.2153    8.2529
        7.6393    9.4764
        1.5041    3.3370
        3.9855    8.3138
        1.8500    5.0079
       1.2631    8.6463
       3.8957    4.9014
       1.9751    8.6199
       1.2565    6.4558
       4.3732    6.1261
       0.4297    8.3551
       3.6931    6.6134
       7.8164    9.8767
       4.8561    8.7376
       6.7750    7.9386
       2.3734    4.7740
       0.8746    3.0892
       2.3088    9.0919
       2.5520    9.0469
       3.3773    6.1886
       0.8690    3.7550
       1.8738    8.1053
       0.9469    5.1476
       0.9718    5.5951
       0.4309    7.5763
       2.2699    7.1371
       4.5000    2.0000  ];

B = [7.2450    3.4422
       7.7030    5.0965
       5.7670    2.8791
       3.6610    1.5002
       9.4633    6.5084
       9.8221    1.9383
       8.2874    4.9380
       5.9078    0.4489
       4.9810    0.5962
       5.1516    0.5319
       8.4363    5.9467
       8.4240    4.9696
       7.6240    1.7988
       3.4473    0.2725
       9.0528    4.7106
       9.1046    3.2798
       6.9110    0.1745
       5.1235    3.3181
       7.5051    3.3392
       6.3283    4.1555
       6.1585    1.5058
       8.3827    7.2617
       5.2841    2.7510
       5.1412    1.9314
       6.0290    1.9818
       5.8863    1.0087
       9.5110    1.3298
       9.3170    1.0890
       6.5170    1.4606
       9.8621    4.3674
       6.0000    8.0000
       2.0000    3.0000  ];

%il training set è l'unione di A e B
T = [A;B];

%dimensione set A e B
nA = size(A,1);
nB = size(B,1);

%plotto i punti
figure(1), clf
scatter(A(:,1),A(:,2), 'MarkerFaceColor','b')
hold on
scatter(B(:,1),B(:,2), 'MarkerFaceColor','r')
title("Points A (blue) and B (red) of the dataset");

%Stavolta la funzione obiettivo è: f = -1/2 l^T X^T X l + e^T l

%{
    Matrice X. Questa è nxL dove n è la dimensione delle componenti di ogni
    punto (2 nel nostro caso), mentre L il numero di punti nel training
    set, ossia nA+nB. Ogni cella (i,j) di X è riempita cosi:
                            y^j * x_i^j
%}

%creo intanto le etichette, 1 per i punti in A, -1 per quelli in B
y = [ ones(nA,1); -ones(nB,1)];

X = zeros(2, nA+nB);

for row=1:size(X,1)
    for column=1:size(X,2)
        X(row,column) = y(column)*T(column,row); 
    end
end

%il vettore e è una serie di 1 
e = ones(1,nA+nB);


% Per utilizzare quadprog dobbiamo trasformarlo in x^T Q x + c^T x
Q = X'*X;

c = e;


%Dx=0, ma abbiamo solo che ly1 + ly2 +... = 0
D = y';

%ogni lambda l deve essere maggiore o uguale a zero
l = zeros(1,nA+nB);

%NUOVO: ogni lambda l deve essere minore o uguale a C
C = 10;
u = C*ones(1,nA+nB);

%{
    risolvo il problema quadratico convesso che siccome era di massimo
    allora inserirò -f
    Nota che Q è stato calcolato non contando -1/2 a sinistra, quindi
    quella parte è come se già avesse contato il "-". Mentre c'è da
    cambiare la parte c^T x.
%}
solution = quadprog(Q,-c,[],[],D,0,l,u,[]);

%la soluzione contiene i vari lambda (pari al numero di punti nel dataset)

%calcolo w1 e w2
w1 = sum(solution.*y.*T(:,1));
w2 = sum(solution.*y.*T(:,2));

%calcolo b
%prendo il primo indice dove lambda è >0 ma (NUOVO) anche <C  <---------- NEW NEW NEW NEW
index = find(solution > 0.001 & solution < C-1e-3);
index = index(1);
%prendo il relativo label e punto
label_index = y(index);
x_index = T(index,:);
b = 1/label_index - [w1 w2]*x_index';

%disegno l'hyperplane w1*x + w2*y + b
max_x = max(T(:,1));

start_x = 0;
start_y = (-w1*start_x-b)/w2;

end_x = max_x;
end_y = (-w1*end_x-b)/w2;

hold on
plot([start_x end_x], [start_y  end_y],'k-')

%evidenzio i punti del dataset che sono vettori di supporto
index = find(solution > 0.001 );
support_vectors = T(index,:);

hold on
scatter(support_vectors(:,1), support_vectors(:,2),100,'d');


