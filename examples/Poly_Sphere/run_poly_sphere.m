%% Create SOS instance
% Number of indeterminates in p
N = 15;
x = msspoly('x', N); % x of dimension N

randn('state',0);

vx = monomials(x,4:4); % Degree 4 homogeneous
% Generate random polynomial
cp = randn(1,length(vx));
cp = (cp > 0) - (cp < 0);
p  = cp*vx;



%% DSOS program

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% DSOS constraint
prog = prog.withDSOS((p - gamma*(x'*x)^2));

% MOSEK options
options = spot_sdp_default_options();
options.solver_options.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solver_options.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.verbose = 0;

% Solve program
tic;
Xddsol = prog.minimize(-gamma, @spot_mosek, options);
dd_toc = toc;

% Optimal value
opt_dsos     = double(Xddsol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_dsos)])

Xdd = double(Xddsol.eval(Xddsol.gramMatrices{1}));





%% SDSOS program

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% DSOS constraint
prog = prog.withSDSOS((p - gamma*(x'*x)^2));

% MOSEK options
options = spot_sdp_default_options();
options.solver_options.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solver_options.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.verbose = 0;

% Solve program
tic;
Xsddsol = prog.minimize(-gamma, @spot_mosek, options);
sdd_toc = toc;

% Optimal value
opt_sdsos     = double(Xsddsol.eval(gamma));

disp(['Optimal value (SDSOS): ' num2str(opt_sdsos)])


%% SOS program

% Using MOSEK for large problems
% is not recommended.
if N <= 25
    
    % Initialize program
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);

    % New free variable gamma
    [prog,gamma] = prog.newFree(1);

    % SOS constraint
    prog = prog.withSOS((p - gamma*(x'*x)^2)); 
    
    % MOSEK options
    options = spot_sdp_default_options();
    
    % Solve program
    tic;
    sol     = prog.minimize(-gamma, @spot_mosek, options);
    sos_toc = toc;
    
    % Optimal value
    opt_sos     = double(sol.eval(gamma));

    disp(['Optimal value (SOS): ' num2str(opt_sos)])

    Xsos = double(sol.eval(sol.gramMatrices{1}));

end



%% Transform problem

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% SOS constraint
prog = prog.withSOS((p - gamma*(x'*x)^2));

% MOSEK options
options = spot_sdp_default_options();

% Build SDP problem data
[pr,~,~,~,~,~] = prog.buildSOSDecompPrimal(prog.sosExpr(1),@newPSD,options);
[P,A,b,c,K,d]  = toSedumi(pr,-gamma);

% Adapt the SeDuMi convention to the SCS convention
A                  = A*P';
c                  = P*c;
A(:,K.f+K.l+1:end) = flip_A(A(:,K.f+K.l+1:end));
c(K.f+K.l+1:end)   = flip_A(c(K.f+K.l+1:end)')';

%% Run Renegar's Method

clear options

if N <= 15
    options.epsilon   = @(i) 1;
    options.useEig      = 1;
else
    options.epsilon   = @(i) 10;
    options.useEig      = 0;
end

options.maxits    = 500;
% options.mu      = 0.02;  % for smoothing
options.printEvery   = 100; % when to print objective
options.smaller_step = 50; % when to decrease step size
options.tol         = 1e-15;
options.epsilon_min = 1e-8; % lower-bound on stepsize
options.cgIter      = 15; % number of conjugate gradient iterations
% options.maxTime = 100;
options.radup   = 0;

ind = triu(ones(K.s)) == 1; % indices for matricising vectors

% Create objective function
%free_map_mat   = pinv( full(A(:,1:K.f)) ) * A(:,K.f+1:end); % map to recover free variables
%free_map       = @(x) -pinv( full(A(:,1:K.f)) ) * (A(:,K.f+1:end) * x - b);
options.objFun = @(x) -c'*x;

% Create initial point satisfying affine constraints
x0         = A\b;
options.x0 = x0;

% if we need strictly feasible point, solve DSOS problem.
% The DSOS solution will not always be full-rank, so check
% failure flag.
%
% obj = p - gamma*(x'*x)^2
% [e,fail_flag] = strictly_feasible_point(prog,obj,A,b);

% Form strictly feasible point.
nk = nchoosek(N+1,2);
C  = ones(nk,1);
j     = 1;
count = 2;
while j <= nk
    C(j) = 1/2;
    j = j + count;
    count = count + 1;
end

C = diag(sparse(C));

% The "50" here is heuristic. Choose a value that makes e positive
% definite. For this problem, choosing different values will not 
% affect the performance of the algorithm (up to the effects of
% numerical conditioning).
e_sdp  = my_smat(x0(K.f+K.l+1:end)) + 50*C;

% create a map to recover free variables from x
free_mat = -pinv(full(A(:,1:K.f+K.l)));
free_map = @(x) free_mat * (A(:,K.f+K.l+1:end) * x - b);

% free variables associated with e_sdp
e_free = free_map(my_svec(e_sdp));
e      = [ e_free; my_svec(e_sdp) ];

[renSol, err] = anySOS(A,b,c,K,e,options);

disp(['Optimal value (AnySOS): ' num2str(err(end,1))])
fprintf('\n')
disp('In Total:')
disp(['Optimal value (AnySOS): ' num2str(err(end,1)) '     Time (AnySOS): ' num2str(err(end,4))])
if N <= 25
    disp(['Optimal value (SOS): ' num2str(opt_sos) '        Time (SOS): ' num2str(sos_toc)])
else
    disp(['Optimal value (SOS): ' 'N/A' '        Time (SOS): ' 'N/A'])
end
disp(['Optimal value (SDSOS): ' num2str(opt_sdsos) '     Time (SDSOS): ' num2str(sdd_toc)])
disp(['Optimal value (DSOS): ' num2str(opt_dsos) '      Time (DSOS): ' num2str(dd_toc)])


