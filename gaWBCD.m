%% Genetic Algorithm

function [a,Fc] = gaWBCD(nbIndiv, nbSelec, nbRules, bp, maxGen)
    tic
    data = csvread('breast-cancer-wisconsin.csv');
    t = 1;

    L = genInPop(nbIndiv, nbRules);
    [f,S] = evalPop(L,nbRules,data);

    while t <= maxGen
        Lsel = selGroup(L,f,S,nbSelec);
        Lchild = crossPop(Lsel,nbIndiv);
        L = mutatePop(Lchild,bp);
        [f,S] = evalPop(L,nbRules,data);
        t = t + 1;
    end

    fmax = f(1); imax = 1;
    for i = 2:nbIndiv
        if f(i) > fmax
            fmax = f(i);
            imax = i;
        end
    end

    a = newFuzzy(L(imax,:),nbRules);
    n = length(data(:,1));
    nbCorrect = 0;
    for i = 1:n
        eval = evalfis(data(i,2:10),a);
        normEv = norm(eval - data(i,11));
        nbCorrect = nbCorrect + (normEv < 1);
    end
    Fc = nbCorrect/n;
    toc
end

%% Create a fuzzy system with the given binary number
function a = newFuzzy(x,Nr)

[P,d,A,w,k] = convert(x,Nr);

a = newfis('WBCD');
a = addvar(a,'input','v1',[1 10]); 
a = addmf(a,'input',1,'low','trapmf',[1 1 P(1) P(1)+d(1)]);
a = addmf(a,'input',1,'high','trapmf',[P(1) P(1)+d(1) P(1)+d(1) P(1)+d(1)]);
a = addvar(a,'input','v2',[1 10]); 
a = addmf(a,'input',2,'low','trapmf',[1 1 P(2) P(2)+d(2)]);
a = addmf(a,'input',2,'high','trapmf',[P(2) P(2)+d(2) P(2)+d(2) P(2)+d(2)]);
a = addvar(a,'input','v3',[1 10]); 
a = addmf(a,'input',3,'low','trapmf',[1 1 P(3) P(3)+d(3)]);
a = addmf(a,'input',3,'high','trapmf',[P(3) P(3)+d(3) P(3)+d(3) P(3)+d(3)]);
a = addvar(a,'input','v4',[1 10]); 
a = addmf(a,'input',4,'low','trapmf',[1 1 P(4) P(4)+d(4)]);
a = addmf(a,'input',4,'high','trapmf',[P(4) P(4)+d(4) P(4)+d(4) P(4)+d(4)]);
a = addvar(a,'input','v5',[1 10]); 
a = addmf(a,'input',5,'low','trapmf',[1 1 P(5) P(5)+d(5)]);
a = addmf(a,'input',5,'high','trapmf',[P(5) P(5)+d(5) P(5)+d(5) P(5)+d(5)]);
a = addvar(a,'input','v6',[1 10]); 
a = addmf(a,'input',6,'low','trapmf',[1 1 P(6) P(6)+d(6)]);
a = addmf(a,'input',6,'high','trapmf',[P(6) P(6)+d(6) P(6)+d(6) P(6)+d(6)]);
a = addvar(a,'input','v7',[1 10]); 
a = addmf(a,'input',7,'low','trapmf',[1 1 P(7) P(7)+d(7)]);
a = addmf(a,'input',7,'high','trapmf',[P(7) P(7)+d(7) P(7)+d(7) P(7)+d(7)]);
a = addvar(a,'input','v8',[1 10]); 
a = addmf(a,'input',8,'low','trapmf',[1 1 P(8) P(8)+d(8)]);
a = addmf(a,'input',8,'high','trapmf',[P(8) P(8)+d(8) P(8)+d(8) P(8)+d(8)]);
a = addvar(a,'input','v9',[1 10]); 
a = addmf(a,'input',9,'low','trapmf',[1 1 P(9) P(9)+d(9)]);
a = addmf(a,'input',9,'high','trapmf',[P(9) P(9)+d(9) P(9)+d(9) P(9)+d(9)]);

a = addvar(a,'output','diagnostic',[2 4]);
a = addmf(a,'output',1,'benign','trapmf',[2 2 2.7 3.3]);
a = addmf(a,'output',1,'malign','trapmf',[2.7 3.3 4 4]);

n = length(A(:,1)); m = length(A(1,:));

for i = 1:n
    for j = 1:m
        if A(i,j) == 3
            A(i,j) = 0;
        end
    end
end

sumAj = sum(A,2); index = [];
for i = 1:n
    if sumAj(i) > 0
        index = [index i];
    end
end

len = length(index);
ruleList = zeros(len+2,m+3);

for i = 1:len
    ruleList(i,1:m) = A(index(i),:);
end

ruleList(1:len,m+1:m+3) = 1;
ruleList(len+1,k) = 1;
ruleList(len+2,k) = 2;
ruleList(len+1:len+2,m+1) = 2;
ruleList(len+1:len+2,m+2) = w;
ruleList(len+1:len+2,m+3) = 1;

a = addrule(a,ruleList);

end

%% Conversion bits into parameters

function [P,d,A,w,k] = convert(x,Nr)
    m = 9;
    P = ones(1,m); A = zeros(Nr,m);
    d = ones(1,m);
    
    for i = 1:1:m
        P(i) = convertBit(x(3*(i-1)+1:3*(i-1)+3));
    end
    
    base = 3*m+1;
    
    for i = 1:1:m
        d(i) = convertBit(x(base+3*(i-1):base+3*(i-1)+2));
    end
    
    base = 6*m+1;
    
    for i = 1:1:Nr
        for j = 1:1:m
            A(i,j) = convertBit(x(base+2*(j-1)+18*(i-1):base+2*(j-1)+18*(i-1)+1));
            if A(i,j) > 2
                A(i,j) = 0;
            end
        end
    end
    
    base = 6*m+2*Nr*m+1;
    w = convertBit(x(base:base+3))/10;
    if w > 1
        w = w - 1;
    end
    
    k = convertBit(x(base+4:base+7));
    if k > 9
        k = k - 9;
    end
end

function y = convertBit(x)
    n = length(x);
    y = 0;
    for i = 0:1:n-1
        y = y + x(n-i)*2^(i);
    end
    if x == zeros(1,n)
        y = 2^n;
    end
end

%% Fitness function

function f = fitness(x, Nr, data)
    alpha = 0.05; beta = 0.01;
    a = newFuzzy(x,Nr);
    
    % Average number of variables per active rule
    nbAnte = 0;
    for i = 1:Nr
        for j = 1:9
            nbAnte = nbAnte + (a.rule(i).antecedent(j) > 0);
        end
    end
    Fv = nbAnte/Nr;
    
    % Classification performance
    n = length(data(:,1));
    nbCorrect = 0; quadDiff = 0;
    for i = 1:n
        eval = evalfis(data(i,2:10),a);
        normEv = norm(eval - data(i,11));
        quadDiff = quadDiff + normEv;
        nbCorrect = nbCorrect + (normEv < 1);
    end
    Fc = nbCorrect/n; Fe = quadDiff/n;
    
    f = Fc - alpha*Fv - beta*Fe;
end

%% Initial population

function L = genInPop(nbIndiv, Nr)
    n = 62+18*Nr;
    L = zeros(nbIndiv, n);
    for i = 1:nbIndiv
        for j = 1:n
            r = rand(1);
            if r > 0.5
                L(i,j) = 1;
            end
        end
    end
end

%% Evaluation of the population
% Return the scaled fitness and the sum of the fitness

function [f,S] = evalPop(L,Nr,data)
    n = length(L(:,1)); f = zeros(n,1);
    imin = 1; imax = 1;
    f(1) = fitness(L(1,:), Nr, data);
    fmin = f(1); fmax = f(1);
    
    for i = 2:1:n
        f(i) = fitness(L(i,:), Nr, data);
        if f(i) < fmin
            fmin = f(i);
            imin = i;
        elseif f(i) > fmax
            fmax = f(i);
            imax = i;
        end         
    end
    
    C = 0.1*fmax - 1.1*fmin;
    D = max(1, fmax + C);
    f = (f + C)/D;
    S = sum(f);
end

%% Select some elements of a population

function Lsel = selGroup(L,f,S,nbSelec)
    Lold = L; Sold = S; fold = f;
    m = length(L(1,:));
    Lsel = zeros(nbSelec,m);
    for i = 1:1:nbSelec
        [x,Lnew,fnew,Snew] = selIndiv(Lold,fold,Sold);
        Lold = Lnew; Sold = Snew; fold = fnew; Lsel(i,:) = x;
    end
end

% Select one element and delete it from the previous population
% Also delete its fitness and decrease the sum of the fitness

function [x,Lnew,fnew,Snew] = selIndiv(L,f,S)
    r = rand(1);
    n = length(L(:,1));
    i = 1; sumFit = f(1);
    while (i < n)&&(sumFit <= r*S)
        i = i + 1;
        sumFit = sumFit + f(i);
    end
    x = L(i,:);
    Snew = S - f(i);
    if i == 1
        Lnew = L(2:n,:);
        fnew = f(2:n);
    elseif i == n
        Lnew = L(1:n-1,:);
        fnew = f(1:n-1);
    else
        L1 = L(1:i-1,:); L2 = L(i+1:n,:); Lnew = [L1; L2];
        f1 = f(1:i-1); f2 = f(i+1:n); fnew = [f1; f2];      
    end
end

%% Create a new population by crossover

function Lchild = crossPop(L,nbIndiv)
    m = length(L(1,:));
    Lchild = zeros(nbIndiv,m);
    for i = 1:floor(nbIndiv/2)
        [Lchild(2*i-1,:),Lchild(2*i,:)] = crossIndiv(L);
    end
    if 2*i < nbIndiv
        [x1,x2] = crossIndiv(L);
        Lchild(nbIndiv,:) = x1;
    end
end

% Create two children by crossover

function [x1,x2] = crossIndiv(L)
    n = length(L(:,1)); m = length(L(1,:));
    a = ceil(n*rand(1)); b = ceil(n*rand(1));
    if a == b
        if b == n
            b = b - 1;
        else
            b = b + 1;
        end
    end
    P1 = L(a,:); P2 = L(b,:); x1 = P1; x2 = P2;
    i = ceil(m*rand(1)); j = ceil(m*rand(1));
    if j > i
        temp = j; j = i; i = temp;
    elseif i == j
        if j == m
            i = i - 1;
        else
            j = j + 1;
        end
    end
    for k = i:1:j
        x1(k) = P2(k); x2(k) = P1(k);
    end
end

%% Mutate a population
function L = mutatePop(L,bp)
    n = length(L(:,1));
    for i = 1:n
        L(i,:) = mutateIndiv(L(i,:),bp);
    end
end

% Mutate an individual
function x = mutateIndiv(x,bp)
    m = length(x);
    for i = 1:m
        r = rand(1);
        if r < bp
            x(i) = mod(x(i)+1,2);
        end
    end
end