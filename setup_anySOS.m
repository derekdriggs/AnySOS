%% Setup Paths
baseDirectory = fileparts(mfilename('fullpath'));
addpath(genpath(baseDirectory));

%% Check Dependencies
if exist('mosek','dir') ~= 7
    warning('MOSEK might not be on path. Demo features require MOSEK.')
end

if exist('spotless-master','dir') ~= 7
    warning('spotless might not be on path. Demo features require spotless.')
end

if exist('cvx','dir') ~= 7
    warning('CVX might not be on path. Some features require CVX.')
end