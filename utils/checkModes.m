function [W_out,Lambda_out,zeroInd] = checkModes(W_in,D_in)

% Algorithm as described in "checkModes" to reorder the modes and return the index of
% the mode with zero eigenvalues
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper, "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E, Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.

%%%%***************************************************************

m = size(W_in, 2);      % number of modes
Lambda_out = D_in;      % initialize the ouput eigenvalues
W_out = W_in;           % initialize the ouput eigenvectors
zeroInd = [];

for k = 1:m
    if D_in(k,k) < 1e-14  % check if the eigenvalue is close to zero
        zeroInd = k;
        for kk = k:m-1
            Lambda_out(kk, kk) = D_in(kk+1, kk+1);
            W_out(:, kk) = W_in(:, kk+1);
        end
        Lambda_out(end, end) = D_in(k, k);
        W_out(:, end) = W_in(:, k);
    end
end
if isempty(zeroInd)
    [diagD,J] = sort(diag(Lambda_out),'descend'); % reorder modes
    W_out = W_out(:, J);
    Lambda_out = diag(diagD);
end