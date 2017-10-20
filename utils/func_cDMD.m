function [Lambda, Phi, W] = func_cDMD(Y1, Y2, C, r, X1, X2, compression, Psi)

% Algorithm as described in "Compressive DMD" to compute compressed DMD &
% compressed sensing DMD
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper: "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E. Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.
% Reference: S. L. Brunton, J. L. Proctor, J. H. Tu, and J. N. Kutz. ?Compressed sensing and
% dynamic mode decomposition?. Journal of Computa- tional Dynamics 2.2 (2015), pp. 165?191.

% Input:
%     measurement matrix Y1, Y2 (shifted),
%     measurement matrix C,
%     sparsifying basis \psi,
%     target rank r,
%     [optional: full state matrix X1, X2].
% Output:
%     cDMD spectrum \LambdaY,
%     cDMD modes \Phi_C or \Phi_CS.


%%%%***************************************************************

% truncated r-rank SVD of Y
[UY, SY, VY] = svd(Y1,'econ');

% low-rank approximation of A_Y
AtildeY = UY(:,1:r)'*Y2*VY(:,1:r)*SY(1:r,1:r)^(-1);

% eigendecomposition of Atilde_Y
[W, Lambda] = eig(AtildeY);

% check zero eigenvalues
[W, Lambda, zeroInd] = checkModes(W,Lambda);

% perform compressed DMD or compressed sensing DMD
switch compression
    case 'C'
        % X is known: perform compressed DMD
        % estimate DMD modes of A
        Phi = X2*VY(:, 1:r)*(SY(1:r, 1:r))^(-1)*W*Lambda^(-1);
        
        % DMD modes for zero eigenvalues
        if ~isempty(zeroInd)
            disp('zero eigenvalue');
            [U, ~, ~] = svd(X1, 'econ');
            Phi(:, end) = U(:, 1:r)*W(:, end);
        end
        
    case 'CS'
        % X is unknown: perform compressed sensing DMD
        % DMD modes of AY
        PhiY = Y2*VY(:, 1:r)*(SY(1:r, 1:r))^(-1)*W*Lambda^(-1);
        
        % DMD modes for zero eigenvalues
        if ~isempty(zeroInd)
            disp('zero eigenvalue');
            PhiY(:, end) = UY(:, 1:r)*W(:, end);
        end
        
        % construct \Theta = C*sparsifying basis
        if nargin > 7
            Theta = C*Psi;
        else
            Theta = dct(C')';
        end
        
        % perform compressed sensing for each mode
        [m, n] = size(C);
        PhiS = zeros(n, r);
        Phi = zeros(n, r);
        for k = 1:r
            % perform l1 minimization on \Phi_Y to solve for \Phi_S
            PhiS(:, k) = cosamp(Theta, PhiY(:, k), 4, m^-6, 10);
            % estimate DMD modes of A
            if nargin > 7 % custom sparsifying basis
                Phi(:, k) = Psi*PhiS(:, k);
            else          % DCT basis
                Phi(:, k) = idct(PhiS(:, k));
            end
        end
end
