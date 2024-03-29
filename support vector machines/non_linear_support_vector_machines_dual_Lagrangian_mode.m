
%{
    Stavolta i dati non sono linearmente separabili.
    
    Risolveremo il duale Lagrangiano per ottenere le soluzioni (lambda) con le quali
    costruire l'hyperplane nel feature space. Useremo come kernel quello gaussiano con 
    gamma=1 per ottenere una superfice più smussata ed evitare l'overfitting.

%}




%{
    DATASET
%}


A = [0.0113    0.2713
        0.9018   -0.1121
        0.2624   -0.2899
        0.3049    0.2100
       -0.2255   -0.7156
       -0.9497   -0.1578
       -0.6318    0.4516
       -0.2593    0.6831
        0.4685    0.1421
       -0.4694    0.8492
       -0.5525   -0.2529
       -0.8250    0.2802
        0.4463   -0.3051
        0.3212   -0.2323
        0.2547   -0.9567
        0.4917    0.6262
       -0.2334    0.2346
        0.1510    0.0601
       -0.4499   -0.5027
       -0.0967   -0.5446  ];

B = [ 1.2178    1.9444
      -1.8800    0.1427
      -1.6517    1.2084
       1.9566   -1.7322
       1.7576   -1.9273
       0.7354    1.1349
       0.1366    1.5414
       1.5960    0.5038
      -1.4485   -1.1288
      -1.2714   -1.8327
      -1.5722    0.4658
       1.7586   -0.5822
      -0.3575    1.9374
       1.7823    0.7066
       1.9532    1.0673
      -1.0233   -0.8180
       1.0021    0.3341
       0.0473   -1.6696
       0.8783    1.9846
      -0.5819    1.8850  ];

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

Q = zeros(nA+nB, nA+nB);


for row=1:size(Q,1)
    for column=1:size(Q,2)
        Q(row,column) = y(row)*y(column)*gaussian_kernel( T(row,:), T(column,:));
    end
end

%il vettore e è una serie di 1 
e = ones(1,nA+nB);


% Abbiamo già la forma per quadprog --> in x^T Q x + c^T x

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


%{
    Non posso disegnarmi l'Hyperplane in quanto questo è lineare solo nel
    feature space. Immagina che lo spazio di partenza sia 2D e che il feature
    space sia 3D. Per trovare la curva (non lineare quindi) nello spazio
    di partenza basta plottare solo markare quei punti (x,y) dove
    l'hyperplane nel feature space è zero, ossia h(x,y)=0.
%}

%prendo il minimo x e y del training set
min_x = min(T(:,1));
min_y = min(T(:,2));
%prendo il massimo di x e y del training set
max_x = max(T(:,1));
max_y = max(T(:,2));
%{
    Per disegnare la curva di boundary, Immagino ogni possibile combinazione 
    tra quei punti creando una sorta di griglia di punti nel rettangolo formato 
    dai punti di min e max sopra e verifico dove h(x)=0
%}

%il termine noto b^* è sempre costante, non dipende da x
%trovo quel lambda tale che il suo valore sia tra 0 e C e prendo il relativo y
index = find(solution > 0.001 & solution < C-1e-3);
index = index(1);
yi = y(index);
xi = T(index,:);
b = 1/yi - sum( solution.*y.*arrayfun(@(row) gaussian_kernel(T(row,:),xi) ,1:(nA+nB) )'  );

%salverò qui i punti g della griglia che fanno da boundary(h(g)=0)
g_boundary = [];
%per ogni punto possibile della griglia...
for g_x=min_x:0.01:max_x
    for g_y=min_y:0.01 :max_y

        %punto della griglia
        g = [g_x,g_y];
        
        %calcolo il w^* che si avrebbe con tale punto
        w = sum(solution.*y.*arrayfun(@(row) gaussian_kernel(T(row,:),g) ,1:(nA+nB) )');

        %calcolo il valore che assume nel feature space
        h_g = w+b;

        %se è zero (con un certo errore) allora posso salvarmelo perchè fa
        %parte del boundary
        if abs(h_g)>1e-2
            g_boundary = [g_boundary; g];
        end
    end
end

%plotto i punti di boundary
hold on
scatter(g_boundary(:,1),g_boundary(:,2), 'MarkerFaceColor','g')



function k = gaussian_kernel(x1, x2)
    gamma = 1;
    k = exp(-gamma*norm(x1-x2)^2);

end



