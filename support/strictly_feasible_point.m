% TODO: output free and nonnegative variables; remove CVX dependence
%
function [e,A,b,c,K] = strictly_feasible_point(A,b,c,K,fr_opts)
% Use CVX to generate a strictly feasible point as a solution to
% a linear program (formed by replacing the SDP with a diagonally-dominant 
% constraint and maximising the smallest diagonal entry). 
% Can also perform facial reduction using fr_lib
% if no strictly feasible point exists.
%
% Inputs:
%   prog - spotsosprog object
%   obj  - msspoly describing SDP objective
%   fr_opts - options for facial reduction using fr_lib TODO: Check if installed
% Outputs:
%   e       - solution to the DSOS problem
%   A,b,c,K - facially reduced problem data (if fr_opts.fr == 1)

if isempty(fr_opts); fr_opts.fr = 0; end

fr    = setOpts(fr_opts, 'fr', 0);

% if performing facial reduction, extract other options
if fr
    cone_type = setOpts(fr_opts, 'cone_type', 'd');
    [A,b,c,K] = facial_reduction(A,b,c,K,cone_type);
end

ind = triu(ones(K.s)) == 1; % indices for matricising vectors

cvx_begin

variable e(K.s,K.s) symmetric
variable xf(K.f,1)
variable xnn(K.l,1)
variable epsilon(1,1)

maximize epsilon
%maximize dot(c(1:K.f),xf) + dot(c(K.f+1:K.f+K.l),xnn) + dot(c(K.f+K.l+1:end),e(ind))

diag(e) - sum(abs(e - diag(diag(e))),2) >= epsilon % diagonally dominant constraint

%e - epsilon*eye(size(e)) >= 0

A*[xf(:);xnn(:);e(ind)] == b

cvx_end

cvx_clear



if any(isnan(e(:))) || any(isinf(e(:)))
    error('CVX could not find a feasible point to the diagonally dominant constrained problem.')
end

if epsilon < 1e-8 && fr_opts.fr == 0
    warning('Could not find a strictly feasible point. Try facial reduction.')
elseif epsilon < 1e-8 && fr_opts.fr == 1
    warning('Could not find a strictly feasible point. Try different facial reduction procedure.')
end


end





function out = setOpts(options, opt, default)
    if isfield(options, opt)
        out = options.(opt);
    else
        out = default;
    end
end % function: setOpts



