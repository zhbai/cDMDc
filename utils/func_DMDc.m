function [Lambda, Phi, Bhat, W, Uhat, Utilde1, Utilde2, Stilde, Vtilde, zeroInd] = func_DMDc(X1, X2, Upsilon, r, rtilde, B)

% Algorithm as described in "DMD with control" to compute DMDc
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper: "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E. Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.
% Reference: J. L. Proctor, S. L. Brunton, and J. N. Kutz. ?Dynamic Mode Decomposition
% with Control?. SIAM J. Applied Dynamical Systems 15.1 (2016), pp. 142?161.

% Input:
%     data matrix X1,
%     shifted data matrix X2,
%     input snapshot matrix \Upsilon,
%     target rank r of X1 or X2 and rtilde of \Omega.

% Output:
%     DMDc spectrum \Lambda,
%     DMDc modes \Phi,
%     [optional: actuation matrix Bhat].

%%%%***************************************************************

% determine whether B is known or unknown
switch nargin
    case 6
        % B is known: perform DMD (adjusted for known actuation).
        [Lambda, Phi, W] = func_DMD(X1, X2-B*Upsilon, r);
        
    case 5
        % B is unknown: perform DMDc
        % matrix of the state and input snapshots
        Omega = [X1; Upsilon];
        
        % truncated rtilde-rank SVD of \Omega
        [Utilde, Stilde, Vtilde] = svd(Omega, 'econ');
        
        % truncated r-rank SVD of X2
        [Uhat, ~, ~] = svd(X2, 'econ');
        
        % split Ustilde into two components
        Utilde1 = Utilde(1:end-rtilde+r, 1:rtilde);
        Utilde2 = Utilde(end-rtilde+r+1:end, 1:rtilde);
        
        % low-rank approximation of A
        temp = Vtilde(:,1:rtilde)*(Stilde(1:rtilde,1:rtilde))^(-1)*Utilde1'*Uhat(:,1:r);
        Atilde = Uhat(:,1:r)'*X2*temp;
        
        % eigendecomposition of Atilde
        [W, Lambda] = eig(Atilde);
        
        % estimate B
        Bhat = X2*Vtilde(:,1:rtilde)*(Stilde(1:rtilde,1:rtilde))^(-1)*Utilde2';
        
        % estimate reduced actuation matrix Btilde
        % Btilde = Uhat'*Bhat;
        
        % check zero eigenvalues
        [W, Lambda, zeroInd] = checkModes(W, Lambda);
        
        % DMD modes of A
        Phi = X2*temp*W;
        
        % DMD modes for zero eigenvalues
        if ~isempty(zeroInd)
            disp('zero eigenvalue');
            Phi(:, end) = Utilde1*Utilde1'*Uhat(:, 1:r)*W(:, end);
        end
end
