
dataset = [
    1.2734    6.2721
2.7453    7.4345
1.6954    8.6408
1.1044    8.6364
4.8187    7.3664
2.7224    6.3303
4.8462    8.4123
4.0497    6.7696
1.0294    8.6174
3.7202    5.1327
3.8238    7.1297
3.5805    7.8660
3.2092    5.7172
1.8724    6.3461
4.0895    5.7509
1.9121    6.2877
2.4835    6.6154
4.5637    7.1943
4.4255    5.1950
2.6097    7.2109
6.0992    8.0496
5.9660    7.3042
5.9726    7.9907
5.6166    7.5821
8.8257    5.4929
8.7426    7.0176
8.2749    6.3890
7.9130    5.3686
5.7032    5.5914
6.4415    5.7927
5.7552    7.6891
5.0048    6.7260
6.2657    7.7776
7.7985    6.0271
7.5010    5.0390
7.1722    7.1291
6.7561    6.1176
6.1497    8.7849
7.0066    8.6258
8.0462    6.5707
3.0994    1.7722
5.6857    2.3666
6.3487    4.7316
6.8860    2.5627
3.2277    2.0929
4.8013    1.6078
5.3299    2.5884
5.7466    2.4989
5.8777    1.5245
5.6002    2.7402
5.9077    1.3661
4.4954    3.4585
5.3263    1.0439
3.4645    3.2930
3.2306    4.1589
6.9191    1.9415
4.1393    2.7921
5.3799    3.2774
6.8486    1.2456
3.7431    2.9852];

%numero di dati
l = size(dataset,1);

%plotto i punti
scatter(dataset(:,1), dataset(:,2), 'black');

%centroidi iniziali
n_centroidi = 3; %quindi 3 cluster
c1 = [5, 7];
c2 = [6, 3];
c3 = [4, 3];
centroids = [c1; c2; c3];

%Abbiamo 3 centroidi e l dati, quindi creo 3*l aij
aij = zeros(l,3);

%per ogni dato pongo aij=1 solo se il centroide j è quello più vicino
for i=1:l
    %prendo il pattern i-esimo
    pattern_i = dataset(i,:);
    %mi salvo la corrente distanza minima
    min_dist = -1;
    %mi salvo l'indice del centroid a cui correntemente è più vicino
    j_near_centroid = -1;

    for j=1:n_centroidi
        
        %calcolo la distanza tra pi e il centroide cj
        centroid = centroids(j,:);

        dist = norm(pattern_i-centroid,1); %<---- Qui la norma è 1

        if min_dist == -1 || min_dist > dist

            min_dist = dist;
            j_near_centroid = j;
        end

    end

    %modifico aij=1 dove i è il pattern i-esimo e j l'indice del centroide più vicino
    aij(i,j_near_centroid)=1;
end

%{
    modifico con colore rosso i punti vicini al centroide c1, verde al
    centroide c2 e blu al centroide c3

%}

%prendo i pattern associati ai centroidi
index_patterns_cluster_1 = find(aij(:,1)==1);
index_patterns_cluster_2 = find(aij(:,2)==1);
index_patterns_cluster_3 = find(aij(:,3)==1);

%li plotto con colori diversi
scatter(dataset(index_patterns_cluster_1,1), dataset(index_patterns_cluster_1,2), 'red');
hold on
scatter(dataset(index_patterns_cluster_2,1), dataset(index_patterns_cluster_2,2), 'green');
scatter(dataset(index_patterns_cluster_3,1), dataset(index_patterns_cluster_3,2), 'blue');

%plotto i centroidi
scatter(c1(1), c1(2),'filled', 'red');
scatter(c2(1), c2(2),'filled', 'green');
scatter(c3(1), c3(2),'filled', 'blue');

hold off
pause(2);
it=1;
current_f_objective = inf;
while true

    fprintf("K-MEANS: Iterazione %d ", it);
    %aggiorno le componenti di ogni cluster come baricentro dei pattern a
    %cui è stato assegnato
    for j=1:n_centroidi
        
        %prendo i pattern a lui assegnati
        index_patterns_cluster_j = find(aij(:,j)==1);
        
        %sommo tutti i pattern tra loro
        numerator = [0.0 0.0];
        for i=1:size(index_patterns_cluster_j,1)

            numerator = numerator + dataset(index_patterns_cluster_j(i),:);
        
        end

        %aggiorno le componenti del centroide corrente come mediana delle
        %componenti dei patter
        %la funzione median ordina già di per se i punti
        %la funzione median calcola la mediana su ogni colonna

        centroids(j,:) = median(dataset(index_patterns_cluster_j,:));

    end

    %Abbiamo 3 centroidi e l dati, quindi creo 3*l aij
    aij = zeros(l,3);

    %per ogni dato pongo aij=1 solo se il centroide j è quello più vicino
    for i=1:l
        %prendo il pattern i-esimo
        pattern_i = dataset(i,:);
        %mi salvo la corrente distanza minima
        min_dist = -1;
        %mi salvo l'indice del centroid a cui correntemente è più vicino
        j_near_centroid = -1;
    
        for j=1:n_centroidi
            
            %calcolo la distanza tra pi e il centroide cj
            centroid = centroids(j,:);
    
            dist = norm(pattern_i-centroid,1);  %<----- Norma 1
    
            if min_dist == -1 || min_dist > dist
    
                min_dist = dist;
                j_near_centroid = j;
            end
    
        end
    
        %modifico aij=1 dove i è il pattern i-esimo e j l'indice del centroide più vicino
        aij(i,j_near_centroid)=1;
    end
    
    %{
        modifico con colore rosso i punti vicini al centroide c1, verde al
        centroide c2 e blu al centroide c3
    
    %}
    
    %prendo i pattern associati ai centroidi
    index_patterns_cluster_1 = find(aij(:,1)==1);
    index_patterns_cluster_2 = find(aij(:,2)==1);
    index_patterns_cluster_3 = find(aij(:,3)==1);
    
    %li plotto con colori diversi
    scatter(dataset(index_patterns_cluster_1,1), dataset(index_patterns_cluster_1,2), 'red');
    hold on
    scatter(dataset(index_patterns_cluster_2,1), dataset(index_patterns_cluster_2,2), 'green');
    scatter(dataset(index_patterns_cluster_3,1), dataset(index_patterns_cluster_3,2), 'blue');
    
    %plotto i centroidi
    scatter(centroids(1,1), centroids(1,2),'filled', 'red');
    scatter(centroids(2,1), centroids(2,2),'filled', 'green');
    scatter(centroids(3,1), centroids(3,2),'filled', 'blue');
    
    %calcolo il valore della funzione obiettivo. Sarebbe la somma delle
    %distanze dei pattern dal proprio centroide
    objective_function_value = 0;
    for i=1:l
        pattern_i = dataset(i,:);
        index_associated_centroid = find(aij(i,:)==1);
        objective_function_value = objective_function_value + norm(dataset(i,:)-centroids(index_associated_centroid,:),1);  %<---- Norma 1
    end

    fprintf("\n  valore funzione obiettivo= %s \n", num2str(objective_function_value));
    title("f(x)= ", num2str(objective_function_value));
    pause(2)
    hold off

    it = it+1;
    if abs(objective_function_value-current_f_objective)<=1e-5
        break
    else
        current_f_objective = objective_function_value;
    end

end

fprintf("fine");
