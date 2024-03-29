
%{

    Questo codice testa se i punti forniti sono punti di minima or weak
    minima in termini di Pareto risolvendo i due rispettivi problemi
    (scalari) ausiliari.

    f(x) = (x1 + 2x2 − 3x3 , −x1 − x2 − x3 , −4x1 − 2x2 + x3)

%}

f1 = @(x) x(1) + 2*x(2) - 3*x(3);
f2 = @(x) -x(1) - x(2) - x(3); 
f3 = @(x) -4*x(1) - 2*x(2) + x(3);

u = [5,0,5];
v = [4,4,2];
w = [1,4,4];

%:::::::::::::::::::: TEST per minimum ::::::::::::::::::::::::::::::::::::

A = [1   2   -3  1  0  0;
     -1 -1   -1  0  1  0;
     -4 -2    1  0  0  1;
     1   1    1  0  0  0;
     0   0    1  0  0  0;];
%per ogni punto da testare calcolo b
bu = [f1(u), f2(u), f3(u),10,5];
bv = [f1(v), f2(v), f3(v),10,5];
bw = [f1(w), f2(w), f3(w),10,5];

%lower bound per epsilon
lb = [0, 0, 0, 0, 0, 0];

f = [0,0,0,1,1,1];

solu = linprog(-f,A,bu,[],[],lb,[]);
solv = linprog(-f,A,bv,[],[],lb,[]);
solw = linprog(-f,A,bw,[],[],lb,[]);

if sum(solu(4:end)==0)
    fprintf("%s è minimum di P\n", mat2str(u));
else
    fprintf("%s NON è minimum di P\n", mat2str(u));
end

if sum(solv(4:end)==0)
    fprintf("%s è minimum di P\n", mat2str(v));
else
    fprintf("%s NON è minimum di P\n", mat2str(v));
end

if sum(solw(4:end)==0)
    fprintf("%s è minimum di P\n", mat2str(w));
else
    fprintf("%s NON è minimum di P\n", mat2str(w));
end


%:::::::::::::::::::: TEST per weak minimum ::::::::::::::::::::::::::::::::::::

A = [1   2   -3  1  0  0 0;
     -1 -1   -1  0  1  0 0;
     -4 -2    1  0  0  1 0;
     1   1    1  0  0  0 0;
     0   0    1  0  0  0 0;
     0   0    0  -1  0  0 1;
     0   0    0   0  -1 0 1;
     0   0    0   0   0 -1 1];

%per ogni punto da testare calcolo b
bu = [f1(u), f2(u), f3(u),10,5,0,0,0];
bv = [f1(v), f2(v), f3(v),10,5,0,0,0];
bw = [f1(w), f2(w), f3(w),10,5,0,0,0];

%lower bound per epsilon
lb = [0, 0, 0, 0, 0, 0,0];

f = [0,0,0,0,0,0,1];

solu = linprog(-f,A,bu,[],[],lb,[]);
solv = linprog(-f,A,bv,[],[],lb,[]);
solw = linprog(-f,A,bw,[],[],lb,[]);

if solu(7)==0
    fprintf("%s è weak minimum di P\n", mat2str(u));
else
    fprintf("%s NON è weak minimum di P\n", mat2str(u));
end

if solv(7)==0
    fprintf("%s è weak minimum di P\n", mat2str(v));
else
    fprintf("%s NON è weak minimum di P\n", mat2str(v));
end

if solw(7)==0
    fprintf("%s è weak minimum di P\n", mat2str(w));
else
    fprintf("%s NON è weak minimum di P\n", mat2str(w));
end