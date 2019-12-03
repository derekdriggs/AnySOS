
function [sol,errHist] = anySOS(A,b,c,K,e,varargin)
% Solves the problem
%
%    min_x   <c,x>   subject to:   A*x = b, x \in C,
%
% where C is the product of the "free" cone, nonnegative cone, or
% positive semidefinite cone.
%
% Inputs:
%   A, b - matrix, vector defining affine contraint
%   c    - matrix in objective function
%   K    - structure containing number of variables in free, nonnegative,
%          and PSD cones in this order.
%   e    - strictly feasible point
%
% Outputs:
%   sol - solutions vector containing free, nonnegative, and semidefinite
%         variables (in that order)
%   errHist - matrix containing objective value, lambda_{min}, lambda_{min}
%             multiple blocks, time, and when the step size is decreased
%
% Options:
%     epsilon    - step size (0.1);
%     maxits     - maximum number of iterations (1000);
%     printEvery - number of iterations between printing (100);
%     objFun     - objective function (c'*x);
%     smaller_step - number of iterations between step size decrease (200);
%     tau          - momentum parameter (0);
%     mu           - smoothing parameter (0);
%     epsilon_min  - minimum step size (1e-12);
%     useEig       - use eig() (use lobpcg() if false) (1);
%     cgIter       - number of conjugate gradient iterations for lobpcg (20);
%     goodEnough   - break after achieving this objective value (-infinity);
%     maxEig       - maximum eigenvalues to use in smoothed gradient (all);
%     clean_up     - make sure solution is feasible (in case of numerical errors) (1);
%     check_feas   - record minimum eigenvalue of each iterate (0);
%     verbose      - print information (1);
%     maxTime      - maximum runtime in seconds (infinity);
%     orthoProj    - use well-conditioned projection onto affine constraints
%                      (use if A is poorly conditioned) (1);
%     radup        - perform radial updates (0);
%     epsilon0     - return to this step size after radial update (1e-1);
%     sol_blocks   - return solution as blocks of free, nonnegative, and
%                      semidefinite varaibles (0);

% Welcome
fprintf([repmat('=',1,76),'\n']);
fprintf('AnySOS by D. Driggs and H. Fawzi -- v1.0\n')
fprintf([repmat('=',1,76),'\n']);

if nargin == 5
    options = [];
elseif nargin == 6
    options = varargin{1};
else
    error('Too many input arguments');
end

% Preprocess the SDP data
[A,b,c,K,e,x0] = preprocess(A,b,c,K,e,options);

nf       = K.f; % number of free variables
nn       = K.l; % number of nonnegative variables
n        = K.s; % number of PSD variables

% Set parameters
epsilon    = setOpts(options, 'epsilon', @(x) 1e-1);
maxits     = setOpts(options, 'maxits', 1000);
printEvery = setOpts(options, 'printEvery',100);
objFun     = setOpts(options, 'objFun', @(x) c'*x );
smaller_step = setOpts(options, 'smaller_step', 200);
tau          = setOpts(options, 'tau',0);
mu           = setOpts(options, 'mu',0);
epsilon_min  = setOpts(options, 'epsilon_min',1e-12);
useEig       = setOpts(options, 'useEig', 1);
cgIter       = setOpts(options, 'cgIter',20);
goodEnough   = setOpts(options, 'goodEnough', []);
maxEig       = setOpts(options, 'maxEig', -1);
clean_up     = setOpts(options, 'clean_up', 1);
check_feas   = setOpts(options, 'check_feas', 0);
verbose      = setOpts(options, 'verbose', 1);
maxTime      = setOpts(options, 'maxTime', []);
orthoProj    = setOpts(options, 'orthoProj',1);
radup        = setOpts(options, 'radup',0);
epsilon0     = setOpts(options, 'epsilon0',@(x) 1e-1);
sol_blocks   = setOpts(options, 'sol_blocks', 0);

% Adjust smoothing parameter
if mu > 0 && mu <= 1e-4
    fprintf('\n Smoothing parameter too small, not using smoothing. \n')
    mu = 0;
end

% Use lobpcg if maxEig is set
if maxEig > 0 && useEig
    fprintf('\n maxEig is set, so using LOBPCG for partial e-val decomposition. \n')
    useEig = 0;
end

% Use eig if maxEigs is too large
if maxEig >= min(n)/5 || n < 25
    fprintf('\n maxEig is large, computing full e-val decomposition. \n')
    useEig = 1;
end

% step size should be function handle
%if isdouble(epsilon)
%    epsilon = @(x) epsilon;
%end

% Additional affine constraint maintains constant objective value
A(end+1,:) = c;

if orthoProj
    [Q,~]      = qr(A',0);
    Qt         = Q';
    
    % Check if Q is sparse (generally Q is not sparse)
    if nnz(Q)/numel(Q) > 0.3
       Q = full(Q);
    end
    
    Proj       = @(x) x - Q*(Qt*x);
    
else
    % Use this method if A is well-conditioned and sparse
    % Create projection onto kernel of A
    P    = inv(A*A');
    Proj = @(x) x - A'*(P*(A*x));
end

% Initialise variables
tic
lambda    = anySOS_line_search(full(x0),full(e),1,K);
x         = e + 1/(1 - lambda)*(x0-e); % project to boundary
x_best    = x;
proj_time = toc;
x_old     = x;        % For momentum
k         = 1;        % Iteration count
obj_last  = c'*x;
obj_best  = obj_last;

% Initialise errHist and record initial objective
errHist      = zeros(maxits,5);
errHist(k,1) = objFun(x); % Starting objective is at e
errHist(k,4) = proj_time;

if verbose
    fprintf('Free variables         : %i                \n',K.f);
    fprintf('Non-negative variables : %i                \n',K.l);
    fprintf('Semidefinite variables : %i block(s) (',max(size(K.s)));
    for i = 1:max(size(K.s))
        fprintf('%ix%i', K.s(i), K.s(i));
        if i ~= max(size(K.s)); fprintf(' '); end
    end
    fprintf(')\n');
    fprintf('Affine constraints     : %i                \n',size(b,1));
    if ~isempty(maxTime)
        fprintf('Stopping after (s)     : %i                \n',maxTime);
    end
    if ~isempty(goodEnough)
        fprintf('Stopping after (obj)   : %i                \n',goodEnough);
    end
    fprintf('Printing every (iter)  : %i                \n',printEvery);
    fprintf([repmat('=',1,76),'\n']);
    fprintf(' iter |   obj   |  time (s) |\n');
    fprintf([repmat('-',1,76),'\n']);
else
    printEvery = -1;
end

tic

while k <= maxits
    x      = x + tau*(x - x_old);
    x_old  = x;
    
    if mu > 0 && maxEig > 0
        % compute gradient of smoothed objective map, using only maxEig
        % eigenvalues (dropping the least important terms).
        [lambda,g] = anySOS_smooth_grad(x,e,mu,K,maxEig,cgIter);
    elseif mu > 0
        [lambda,g] = anySOS_smooth_grad(x,e,mu,K); % compute gradient of smoothed objective
    else
        [lambda,g] = anySOS_subgrad(x,e,K,useEig,maxEig,cgIter);
    end % if: subgradient computation
    
    pik  = e + 1/(1 - lambda) * (x-e);  % project to boundary
    obj  = c'*pik; % compute objective for various tasks
    
    if obj <= obj_best  % store the best iterate
        obj_best = obj;
        x_best = pik;
        if ~isempty(goodEnough) && obj <= goodEnough % break condition
             break
        end % if: good enough solution
        if ~isempty(maxTime) && errHist(k,4) > maxTime % break condition
             break
        end % if: time limit exceeded
    end % if: store best iterate
    
    if printEvery > -1 && mod(k,printEvery) == 0
        fprintf('%5d | %7.2e | %7.2e |\n',...
            k,objFun(pik),errHist(k,4))
    end % if: print
    
    % Decrease the step size if we are not making progress
    if smaller_step > -1 && mod(k,smaller_step) == 0
        if epsilon(k) >= epsilon_min
            epsilon = @(x) epsilon(x)/2; % halve the step size
            if ~useEig
                cgIter = cgIter + 2; % more accurate subgradients
            end % if: cgIter update
            errHist(k+1,5) = 1; % record when we halve the step size
        end % if: objective progress
        % obj_last = obj;
    end % if: smaller step
    
    % Perform a radial update if conditions are met
     if radup && c'*e - obj >= 4/3*c'*(e - x)
         x     = pik;
         x_old = x;    % Reset momentum
         epsilon = @(x) epsilon0(x); % Reset step size 
     else
        g     = Proj( g ); % project onto kernel of A
        g     = g/norm(g(:)); % normalize
        
        alpha = epsilon(k); % step size

        x  = x  + alpha*g;  % subgradient update
    end % if: radial update
    
    errHist(k+1,4) = toc + errHist(k,4);
    errHist(k+1,2) = lambda;
    errHist(k+1,1) = objFun(pik);
    
    if check_feas
        nvec = (n + 1) .* (n ./ 2);
        eig_list = zeros(size(n,2),1);
        for i = 1:size(n,2)
            pik_full    = my_smat(pik(nf+nn+sum(nvec(1:i-1))+1:nf+nn+sum(nvec(1:i))));
            eig_list(i) = min(eig(pik_full));
        end
        errHist(k+1,3) = min(eig_list);
    end
    
    tic
    
    k = k + 1;
    
end % while: main loop

% perform full eigendecomposition for post-processing
if clean_up
    lambda      = anySOS_line_search(full(x_best),full(e),1,K); 
    x_best      = e + 1/(1 - lambda) * (x_best-e);
end % if: clean-up

% If sol_blocks, form solution into blocks of free variables,
%   non-negative variables, and semidefinite variables.
if sol_blocks
    sol.free_var = x_best(1:nf);
    sol.nn_var   = x_best(nf+1:nf+nn);

    nvec = (n + 1) .* (n ./ 2);
    for i = 1:size(n,2)
        xs{i}  = my_smat( x_best(nf+nn+sum(nvec(1:i-1))+1:nf+nn+sum(nvec(1:i))) );
    end

    sol.xs = xs;
else
    sol = x_best;
end

errHist = errHist(1:k,:);

fprintf([repmat('=',1,76),'\n']);

end % function: anySOS




function out = setOpts(options, opt, default)
    if isfield(options, opt)
        out = options.(opt);
    else
        out = default;
    end
end % function: setOpts

