
%{

    Questo script serve solo a testare se un punto è un Nash equilibrium
    sfruttando la proprietà che la gap function faccia zero quando tale
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

%{
    la gap function è fatta dai seguenti termini:
        x^T (C1+C2)y - min (u^T C_1 y) - min (v^T C2 x)
    
    Il primo è una costante, mentre gli altri due sono problemi lineari da
    risolvere rispetto u (il primo) e v (il secondo)
    
%}

constant = x'*(C1+C2)*y;

%:::::::::::::::::::::::: primo problema ::::::::::::::::::::

% u1+u2+u3 = 1
Aeq = [1 1 1];
beq = [1];

%ui>=0
lb = [0 0 0];

c = C1*y;

[sol1,val1] = linprog(c,[],[],Aeq,beq,lb,[]);


%:::::::::::::::::::::::: secondo problema ::::::::::::::::::::

% u1+u2+u3 = 1
Aeq = [1 1];
beq = [1];

%ui>=0
lb = [0 0];

c = x'*C2;

[sol2,val2] = linprog(c,[],[],Aeq,beq,lb,[]);



%testiamo se è un Nash equilibrium
if (constant - val1 - val2)>1e-3
    fprintf("la gap function ha valore %s. Quindi non è un Nash equilibrium.\n",num2str(constant - val1 - val2));
else
    fprintf("la gap function ha valore %s. Quindi è un Nash equilibrium.\n",num2str(constant - val1 - val2));
end