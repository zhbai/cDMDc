function [Lambda, Phi, W] = func_DMD(X1, X2, r)

% Algorithm as described in "Exact DMD" to compute DMD
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper: "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E. Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.
% Reference: J. H. Tu, C. W. Rowley, D. M. Luchtenburg, S. L. Brunton, and J. N. Kutz. ?On dynamic mode
% decomposition: Theory and applications?. Journal of Computational Dynamics 1.2 (2014), pp. 391?421.

% Input:
%     data matrix X1,
%     shifted data matrix X2,
%     target rank r.
% Output:
%     DMD spectrum \Lambda,
%     DMD modes \Phi.

%%%%***************************************************************

% truncated r-rank SVD of X
[U, S, V] = svd(X1,'econ');
Sinv = S(1:r, 1:r)^(-1);

% low-rank approximation of A
Atilde = U(:, 1:r)'*X2*V(:, 1:r)*Sinv;

% eigen-decomposition of Atilde
[W,Lambda] = eig(Atilde);

% check zero eigenvalues
[W, Lambda, zeroInd] = checkModes(W, Lambda);

% DMD modes of A
Phi = X2*V(:, 1:r)*Sinv*W*Lambda^(-1);

% DMD modes for zero eigenvalues
if ~isempty(zeroInd)
    disp('zero eigenvalue');
    Phi(:, end) = U(:, 1:r)*W(:, end);
end
