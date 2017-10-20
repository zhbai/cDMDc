function [a,b] = tune_sign(a,b,type)

% Algorithm as described in "tune_sign" to regulate the DMD modes for
% plots and error report.
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper, "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E, Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.

%%%%***************************************************************

m = size(a,2);          % number of modes

if type == 1            % check the sign of the first element of two vectors
    for k = 1:m

        if sign(real(a(1,k))*real(b(1,k))) == -1
            b(:,k) = -b(:,k);
        end
        
    end
elseif type == 2        % check the inner product of the two vectros
    for k = 1:m

        if real(a(:,k))'*real(b(:,k)) < 0
            b(:,k) = -b(:,k);
        end
        
    end
end