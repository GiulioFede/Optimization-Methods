
%{

    Il problema convesso è il seguente:
        min (x1 , x1^2 + x2^2 − 2x1)
        -x1<=0
        -x2<=0
        x1+x2<=2

%}

a1 = 0;
a2 = 0;

A = [-1 0;
      0 -1;
      1  1];

b = [0 0 2];

figure(2), clf
title("Weak Minimum of P");
%:::::::::::::::::::::::::: WEAK MINIMUM ::::::::::::::::::::::::::::::::::
for a1=0:0.01:1

    a2 = 1-a1;

    H = [a1  0;
     0   a2];
    
    H = H*2; %perchè nella definizone c'è 1/2 --> 1/2 x^T H x

    c = [a1-2*a2,0];

    sol = quadprog(H,c,A,b,[],[]);

    fprintf("Soluzione con a1=%s e a2=%s pari a: %s\n", num2str(a1), num2str(a2), mat2str(sol));
    plot(sol(1),sol(2),'g.','MarkerSize',20);
    hold on

end

plot([0 0 2,0], [0 2 0,0]);


figure(3), clf
title("Minimum of P");
%:::::::::::::::::::::::::: WEAK MINIMUM ::::::::::::::::::::::::::::::::::
for a1=0.001:0.001:0.999

    a2 = 1-a1;

    H = [a1  0;
     0   a2];
    
    H = H*2; %perchè nella definizone c'è 1/2 --> 1/2 x^T H x

    c = [a1-2*a2,0];

    sol = quadprog(H,c,A,b,[],[]);

    fprintf("Soluzione con a1=%s e a2=%s pari a: %s\n", num2str(a1), num2str(a2), mat2str(sol));
    plot(sol(1),sol(2),'g.','MarkerSize',20);
    hold on

end

plot([0 0 2,0], [0 2 0,0]);


