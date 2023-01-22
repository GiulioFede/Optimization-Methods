%{

    Questo script serve solo a testare se un punto è un Nash equilibrium
    sfruttando la proprietà che la D gap function faccia zero quando tale
    punto (x,y) dato la rende zero.

%}

%punto da testare se è o meno un Nash Equilibrium

x = [1/3; 1/3; 1/3];
y = [1/2; 1/2];

C1 = [3 3;
      4 1;
      6 0];

C2 = [3 4;
      4 0;
      3 5];

a = 1;
b = 10;

%{
    la gap function è fatta dai seguenti termini:
        x^T (C1+C2)y - a/2*|x|^2 - a/2*|y|^2 - min ( a/2|u|^2 + u*(C1y-ax))
        - min (a/2|v| + v(C2x-ay))

    I primi 3 sono costanti, mentre gli altri due sono problemi quadratici da
    risolvere rispetto u (il primo) e v (il secondo)
    
    La D-gap è invece:
            gap_a(x,y) -gap_b(x,y)

    
%}

constant1 = x'*(C1+C2)*y - a/2*norm(x)^2 -a/2*norm(y)^2;


%:::::::::::::::::: Primo problema quadratico gap_a :::::::::::::::::::::::::::::

H = a*eye(3);
c = C1*y-a*x;

% u1+u2+u3 = 1
Aeq = [1 1 1];
beq = [1];

%ui>=0
lb = [0 0 0];

[sol1,val1] = quadprog(H,c,[],[],Aeq,beq,lb);


%:::::::::::::::::: Secondo problema quadratico gap_a :::::::::::::::::::::::::::::

H = a*eye(2);
c = C2'*x-a*y;

% v1+v2 = 1
Aeq = [1 1];
beq = [1];

%ui>=0
lb = [0 0 0];

[sol3,val2] = quadprog(H,c,[],[],Aeq,beq,lb);





constant2 = x'*(C1+C2)*y - b/2*norm(x)^2 -b/2*norm(y)^2;
%:::::::::::::::::: Primo problema quadratico gap_b :::::::::::::::::::::::::::::

H = b*eye(3);
c = C1*y-b*x;

% u1+u2+u3 = 1
Aeq = [1 1 1];
beq = [1];

%ui>=0
lb = [0 0 0];

[sol4,val3] = quadprog(H,c,[],[],Aeq,beq,lb);


%:::::::::::::::::: Secondo problema quadratico gap_b :::::::::::::::::::::::::::::

H = b*eye(2);
c = C2'*x-b*y;

% v1+v2 = 1
Aeq = [1 1];
beq = [1];

%ui>=0
lb = [0 0 0];

[sol2,val4] = quadprog(H,c,[],[],Aeq,beq,lb);


%valore della D-gap

d_gap_value = (constant1 - val1 - val2) - (constant2 - val3 - val4);



%testiamo se è un Nash equilibrium
if d_gap_value>1e-3
    fprintf("la d-gap function ha valore %s. Quindi non è un Nash equilibrium.\n",num2str(d_gap_value));
else
    fprintf("la d-gap function ha valore %s. Quindi è un Nash equilibrium.\n",num2str(d_gap_value));
end