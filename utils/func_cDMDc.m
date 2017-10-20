function [Lambda, Phi, Bhat] = func_cDMDc(Y1, Y2, C, Upsilon, r, rtilde, X1, X2, B, type, Psi)

% Algorithm as described in "Compressive DMDc" to compute compressed DMDc &
% compressed sensing DMDc
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper, "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E, Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.

% Input:
%     measurement matrix Y1, Y2 (shifted),
%     measurement matrix C,
%     sparsifying basis \psi,
%     target rank r,
%     [optional: full state matrix X1, X2, actuation matrix B].
% Output:
%     cDMD spectrum \LambdaY,
%     cDMD modes \Phi_C or \Phi_CS.

% Four types explored depending on the pre-knowledge of X and B:
% type 1: X is known, B is known;
% type 2: X is knwon, B is unknown;
% type 3: X is unknown, B is known;
% type 4: X is unknown, B is unknown.

%%%%***************************************************************
switch type
    case '1'
        % X is known, B is known: perform compressed DMD
        disp('Case 1: X is known, B is known')
        [Lambda, Phi] = func_cDMD(Y1, Y2 - C*B*Upsilon, C, r, X1, X2-B*Upsilon, 'C');
        
    case '2'
        % X is knwon, B is unknown: perfom DMDc
        disp('Case 2: X is known, B is unknown')
        [Lambda, ~, ~, W, UhatY, Utilde1Y, Utilde2Y, StildeY, VtildeY, zeroIndY] = func_DMDc(Y1, Y2, Upsilon, r, rtilde);
        
        % estimate DMD modes of A
        Phi = X2*VtildeY(:, 1:rtilde)*(StildeY(1:rtilde, 1:rtilde))^(-1)*Utilde1Y'*UhatY(:, 1:r)*W*Lambda^(-1);
        
        % Estimate actuation matrix Bhat
        Bhat = X2*VtildeY(:, 1:rtilde)*(StildeY(1:rtilde, 1:rtilde))^(-1)*Utilde2Y';
        
        % DMD modes for zero eigenvalues
        Omega = [X1; Upsilon];
        [Utilde, ~, ~] = svd(Omega, 'econ');
        Utilde1 = Utilde(1:end-rtilde+r,  1:rtilde);
        if ~isempty(zeroIndY)
            Phi(:,end) = Utilde1*Utilde1Y'*UhatY(:, 1:r)*W(:, end);  % zero eigenvalue
        end
        
    case '3'
        % X is unknown, B is known: perform compressed sensing DMD
        disp('Case 3: X is unknown, B is known')
        [Lambda, Phi] = func_cDMD(Y1, Y2 - C*B*Upsilon, C, r, [], [], 'CS');
        
    case '4'
        % X is unknown, B is unknown: perform DMDc and then reconstruct
        % modes
        disp('Case 4: X is unknown, B is unknown')
        [Lambda, PhiY, BYhat] = func_DMDc(Y1, Y2, Upsilon, r, rtilde);
        [m, n] = size(C);
        PhiS = zeros(n, r);
        Phi = zeros(n, r);
        
        % construct \Theta = C*sparsifying basis
        if nargin > 10      % custom sparsifying basis
            Theta = C*Psi;
        else                % DCT basis
            Theta = dct(C')';
        end
        
        % perform compressed sensing for each mode
        for k = 1:r
            PhiS(:, k) = cosamp(Theta, PhiY(:, k), 4, m^-6, 10);
            if nargin > 10      % custom sparsifying basis
                Phi(:,k) = Psi*PhiS(:,k);
            else                % DCT basis
                Phi(:,k) = idct(PhiS(:,k));
            end
        end
        
        % perform compressed sensing for B
        for k = 1:size(BYhat,2)
            BS = cosamp(Theta, BYhat(:, k), 4, m^-6, 10);
            if nargin > 10      % custom sparsifying basis
                Bhat = Psi*BS(:,k);
            else                % DCT basis
                Bhat = idct(BS(:,k));
            end
        end
end