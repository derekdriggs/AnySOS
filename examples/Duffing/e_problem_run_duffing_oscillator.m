%% Create Network of Duffing Oscillators
% Number of indeterminates in network
N = 10;
x = msspoly('x', N); % x of dimension N

randn('state',0);

p = 0;

bduff = zeros(N);    % initialise coefficients in network
aduff = zeros(N,1);

% Create Hamiltonian function for Duffing oscillator network
for k = 1:N
    
    aduff(k) = rand(1)*(1.5 - 0.5) + 0.5;
    p = p + aduff(k)*(1/2*x(k)^2 - 1/4*x(k)^4);
    
    for i = 1:N
        bduff(i,k) = rand(1)*(1.5/N - 0.5/N) + 0.5/N;
        p = p + bduff(i,k)/8*(x(i) - x(k))^4;
    end
    
end


%% SOS program

if N <= 20 % Be careful with large N

% MOSEK options
options = spot_sdp_default_options();
options.solver_options.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solver_options.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.verbose = 0;
options.scale_monomials = false;

a   = 0;
b   = 2;

tic
for i = 1:6 % bisection loop
    
    g = (a + b)/2;
    
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,lambda] = prog.newFree(N,1);
    [prog,mu] = prog.newFree(1);

    const = lambda'*x.^2*g - lambda'*x.^4 - mu*ones(1,N)*x.^4;

    % SOS constraint
    prog = prog.withSOS(p - const);

    % Solve program
    sol = prog.minimize(mu, @spot_mosek, options);

    feas = sign(double(sol.eval(mu)));
    if feas > 0
        b = g;
    else
        a = g;
    end
        
end

fprintf('SOS value: %6.2f \n', a)

end






%% SDSOS program

% MOSEK options
options = spot_sdp_default_options();
options.solver_options.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solver_options.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.verbose = 0;
options.scale_monomials = false;

a   = 0;
b   = 2;

tic
for i = 1:6
    
    g = (a + b)/2;
    
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,lambda] = prog.newFree(N,1);
    [prog,mu] = prog.newFree(1);

    const = lambda'*x.^2*g - lambda'*x.^4 - mu*ones(1,N)*x.^4;

    % SDSOS constraint
    prog = prog.withSDSOS(p - const);

    % Solve program
    sol = prog.minimize(mu, @spot_mosek, options);
    
    feas = sign(double(sol.eval(mu)));
    if feas > 0
        b = g;
    else
        a = g;
    end
    
end

fprintf('SDSOS value: %6.2f \n', a)



%% DSOS program

% MOSEK options
options = spot_sdp_default_options();
options.solver_options.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solver_options.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.verbose = 0;
options.scale_monomials = false;

a   = 0;
b   = 2;
tic
for i = 1:6
    
    g = (a + b)/2;
    
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,lambda] = prog.newFree(N,1);
    [prog,mu] = prog.newFree(1);

    const = lambda'*x.^2*g - lambda'*x.^4 - mu*ones(1,N)*x.^4;

    % DSOS constraint
    prog = prog.withDSOS(p - const);
    
    % Solve program
    sol = prog.minimize(mu, @spot_mosek, options);
 
    feas = sign(double(sol.eval(mu)));
    if feas > 0
        b = g;
    else
        a = g;
    end
    
end

fprintf('DSOS value: %6.2f \n', a)


%% Form strictly feasible point for AnySOS

% Could try DSOS solution as strictly feasible point
%[isit,e,~] = isDSOS(p - ( zeros(1,N)*x.^2*g - zeros(1,N)*x.^4 - 2*ones(1,N)*x.^4 ) );
%e          = double(e);

tic

% Form strictly feasible point manually using basis in P
prog = spotsosprogtransform;
prog = prog.newIndeterminate('x',N);
[prog,lambda] = prog.newFree(N,1);
[prog,mu] = prog.newFree(1);

prog  = prog.withSOS(p);

% Construct basis
[~,~,~,~,~,P] = prog.transform(mu, @spot_mosek, []);

P = P{1};

M = zeros(nchoosek(N+2,2) - 1);
Pos = M;

for i = 1:N
    
    e = ones(nchoosek(N+2,2) - 1,1);
    ind = find(P-x(i));
    e(ind) = 0;
    E1(:,i) = e;
    
    M = M+aduff(i)/2*E1(:,i)*E1(:,i)';
    
    e = ones(nchoosek(N+2,2) - 1,1);
    ind = find(P-x(i)^2);
    e(ind) = 0;
    
    E2(:,i) = e;
    
    M = M - aduff(i)/4*E2(:,i)*E2(:,i)';
    
    Pos = Pos + 1/4*E2(:,i)*E2(:,i)';
end




for i = 1:N
    for j = 1:N
        e = ones(nchoosek(N+2,2) - 1,1);
        ind = find(P-x(i)*x(j));
        e(ind) = 0;

        E12{i}(:,j) = e;
        
        bvec = E2(:,i) + E2(:,j) - 2*E12{i}(:,j);
        
        M = M+1/8*bduff(i,j)*bvec*bvec';
        
    end

end

e = M + 2*Pos;

clear E2 E1 E12 ind

time = toc;


%% Run Renegar's Method

if N >= 25
    params.maxits    = 10000;
    params.mu        = 0.2;  % for smoothing
    params.epsilon_min = 0.005;
    params.smaller_step = 10; % when to decrease step size
else
    params.maxits    = 10000;
    params.mu        = 0;  % for smoothing
    params.epsilon_min = 0.001;
    params.smaller_step = 50; % when to decrease step size
end

params.printEvery  = 50; % when to print objective
params.tol         = 1e-15;
params.goodEnough  = 0;  % break if objective is below this value

if N >= 25
    params.useEig       = 0;
    params.maxEig       = 5;
    params.cgIter       = 5; % number of conjugate gradient iterations
else
    params.useEig       = 1;
end

params.epsilon = @(x) 0.005; % step size
a   = 0;
b   = 2;


for i = 1:6
    
    tic
    
    g = (a + b)/2
    
    prog = spotsosprogtransform;
    prog = prog.newIndeterminate('x',N);
    [prog,lambda2] = prog.newFree(N,1);
    [prog,mu] = prog.newFree(1);

    const = lambda2'*x.^2*g - lambda2'*x.^4 - mu*ones(1,N)*x.^4;
    prog  = prog.withSOS(p - const);
    
    % Transform problem into SeDuMi format
    [A,y,c,K,~] = prog.transform(mu, @spot_mosek, options);
    
    c(K.f+1:end) = flip_A(c(K.f+1:end)')';
    A(:,K.f+1:end) = flip_A(A(:,K.f+1:end));
    
    ind = triu(ones(K.s)) == 1; % indices for matricising vectors
    
    if i >= 2
        j = 1;
        k = 1;
        inds = zeros(size(P,1),1);
        lambda = flip(xf(1:end-1));

        while j <= size(P,1)
            inds(j) = lambda(k);
            k = k+1;
            j = j+k;
        end
    
        params.x0 = renSol - (g - g_old)*diag(inds);
    else
        % Create initial point satisfying affine constraints
        params.x0 = M - Pos;
        
        clear M
    end

    time = time + toc;
    
    % Call anySOS
    [Rsol, errHist] = anySOS(A,y,c,K,my_svec(e),params);
    xf     = Rsol.free_var;
    renSol = Rsol.xs{1};
    
    time = time + max(errHist(:,4));
    
    feas = sign(xf(end));
    
    if feas > 0
        b = g;
    else
        a = g;
    end
    
    params.mu          = params.mu/2;
    num_smaller        = nnz(errHist(:,end));
    new_eps            = params.epsilon(1)/(2^num_smaller);
    params.epsilon     = @(x) new_eps;
    g_old = g;
    
end



fprintf('Renegar value: %6.2f \n', a)


