function outPhi = normalize(inPhi)

% Algorithm as described in "normalize" to normalize the DMD modes for
% calibration.
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper, "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E, Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.

%%%%***************************************************************

m = size(inPhi,2);      % number of modes

for k = 1:m
    inPhi(:,k) = inPhi(:,k)/inPhi(1,k);  % normalize the first element to 1
    outPhi(:,k) = inPhi(:,k)/norm(inPhi(:,k)); 
end

