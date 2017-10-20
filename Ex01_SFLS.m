% Algorithm as described in "Compressive DMDc" to compute compressed DMDc &
% compressed sensing DMDc
% Copyright 2017, All Rights Reserved
% Code by Zhe Bai
% For Paper, "Dynamic mode decomposition for compressive system identification"
% by Z. Bai, E. Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.

% Example: stochastically forced linear system

% Actuation matrix B:
% 1. in the span of P (projection matrix);
% 2. random generated;
% 3. not in the subspace of P.

% Measurement matrix C: 
% 1. uniform distribution; 
% 2. gaussian distribution;
% 3. single pixel measurement

% Four cases tested:
% 1. X known, B known
% 2. X known, B unknown
% 3. X unknown, B known
% 4. X unknown, B unknown

%%%%***************************************************************

clear all, close all, clc
addpath('./utils');
figpath = './figures/';
outpath = './output/';

%% generate data
n = 1024;               % dimension of states
p = 128;                % dimension of compressed measurements
dt = 0.1;
tspan = [0:dt:30];      % time span
nt = numel(tspan);      % number of time steps
rng(1);                 % fix random generator
Atilde = [0.9  0.2; -0.1  0.9]; % dynamics
Btilde = [0.1; 0.01];   % actuation matrix

for B_choice = 1     % 1-sub, 2-randn or 3-nonsub
    
    for CType = 1    % 1-unifrom, 2-gaussian 3-single pixel
    % get two-dimensional system
    xtilde = [0.25; 0.25];   % initial condition
    Upsilon = zeros(1, nt);
    for tk = 1:nt
        Upsilon(tk) = randn; % generate random input vector
        xtilde(:, tk+1) = Atilde*xtilde(:,tk) + Btilde*(Upsilon(tk));
    end
    
    % projection matrix
    p1 = zeros(n, 1); p1(3) = 1; p1(29) = 1; p1 = p1/norm(p1);
    p2 = zeros(n, 1); p2(11) = 1.5; p2(47) = 1; p2 = p2/norm(p2);
    P = idct2([p1 p2]); % orthogonal columns of P
    
    % true A
    A = P*Atilde*(pinv(P));
    
    % generate B matrix: three cases
    if B_choice == 1
        B = P*Btilde; % 1. span of P
        B_name = '_sub';
    end
    if B_choice == 2
        B = randn(n,1); B = B/norm(B); % 2. random
        B_name = '_randn';
    end
    
    p3 = zeros(n,1); p3(19) = 1;
    if B_choice == 3
        B = idct2(p3); B = B/norm(B); % 4. not span of P
        B_name = '_nonsub';
    end
    
    %% project to get high-dimensional system
    X = P*xtilde(:,1); % initialize X
    for tk = 1:nt
        X(:,tk+1) = A*X(:,tk) + B*Upsilon(tk);
    end
    
    %% real dynamics
    [T,DA] = eig(A);
    T = T(:, 1:2);
    DA = DA(1:2, 1:2);
    
    %% DMDc
    if B_choice == 1  % B is in the span of P
        r = 2;
    else
        r = 3;        % B is not in the span of P
    end
    rtilde = r + 1;
    X1 = X(:, 1:end-1);
    X2 = X(:, 2:end);
    
    % B known
    [D01, Phi01] = func_DMDc(X1, X2, Upsilon, r, rtilde, B);
    
    % B unknown
    [D02, Phi02, Bhat02] = func_DMDc(X1, X2, Upsilon, r, rtilde);
    
    %% cDMDc
    % unifrom distribution
    if CType == 1
        C = randn(p,n);
        % grassian distribution
    elseif CType == 2 
        C = rand(p,n);
        %   single pixel measurement
    elseif CType == 3
        C = zeros(p,n);
        ind = randperm(n);
        for jj = 1:p
            C(jj, ind(jj)) = 1;
        end
    end
    
    %% compressed measurements
    Y = C*X; % compress X to get Y
    Y1 = Y(:, 1:end-1);
    Y2 = Y(:, 2:end); % shift matrix
    
    %% test four cases
    % 1. X known, B known
    [D11, Phi11] = func_cDMDc(Y1, Y2, C, Upsilon, r, rtilde, X1, X2, B, '1');
    
    % 2. X known, B unknown
    [D12, Phi12, Bhat12] = func_cDMDc(Y1, Y2, C, Upsilon, r, rtilde, X1, X2, [], '2');
    
    % 3. X unknown, B known
    [D21, Phi21] = func_cDMDc(Y1, Y2, C, Upsilon, r, rtilde, X1, X2, B, '3');
    
    % 4. X unknown, B unknown
    [D22, Phi22, Bhat22] = func_cDMDc(Y1, Y2, C, Upsilon, r, rtilde, [], [], [], '4');
    
    %% normalize modes
    T = normalize(T);
    Phi01 = normalize(Phi01); Phi02 = normalize(Phi02);
    Phi11 = normalize(Phi11); Phi12 = normalize(Phi12);
    Phi21 = normalize(Phi21); Phi22 = normalize(Phi22);
    [~, Phi01] = tune_sign(T, Phi01, 1); [~, Phi02] = tune_sign(T, Phi02, 1);
    [~, Phi11] = tune_sign(T, Phi11, 1); [~, Phi12] = tune_sign(T, Phi12, 1);
    [~, Phi21] = tune_sign(T, Phi21, 1); [~, Phi22] = tune_sign(T, Phi22, 1);
    if r == 3
        [~, Phi11(:,3)] = tune_sign(Phi01(:,3), Phi11(:,3), 1); [~, Phi12(:,3)] = tune_sign(Phi02(:,3), Phi12(:,3), 1);
        [~, Phi21(:,3)] = tune_sign(Phi01(:,3), Phi21(:,3), 1); [~, Phi22(:,3)] = tune_sign(Phi02(:,3), Phi22(:,3), 1);
    end
    err_Phi_01 = norm(Phi01(:,1:2)-T,'fro')/norm(T,'fro'); err_Phi_02 = norm(Phi02(:,1:2)-T,'fro')/norm(T,'fro');
    err_Phi_11 = norm(Phi11(:,1:2)-T,'fro')/norm(T,'fro'); err_Phi_12 = norm(Phi12(:,1:2)-T,'fro')/norm(T,'fro');
    err_Phi_21 = norm(Phi21(:,1:2)-T,'fro')/norm(T,'fro'); err_Phi_22 = norm(Phi22(:,1:2)-T,'fro')/norm(T,'fro');
    %% plot DMD modes
    for ii = 1:2
        B_cond = num2str(ii);
        for k = 1:2
            figure
            if k == 2 && real(D01(1,1)) == real(D01(2,2))
                plot(imag(T(:,k-1)),'color',[0.5 0.5 0.5], 'linewidth',8);
                hold on
                plot(imag(eval(strcat('Phi0', B_cond,'(:,k-1)'))), 'k-', 'linewidth',4);
                plot(imag(eval(strcat('Phi1', B_cond,'(:,k-1)'))), 'b--', 'linewidth',2);
                plot(1:8:n,imag(eval(strcat('Phi2', B_cond,'(1:8:n,k-1)'))), 'ro-', 'markersize', 4, 'linewidth', 1.2);
                set(gca,'Xtick',[0,256,512,768,1024])
            else
                if k < 3
                    plot(real(T(:,k)), 'color', [0.5 0.5 0.5], 'linewidth',8);
                    hold on,
                end
                plot(real(eval(strcat('Phi0', B_cond,'(:,k)'))),'k-','linewidth',4);
                hold on
                plot(real(eval(strcat('Phi1', B_cond,'(:,k)'))), 'b--','linewidth',2);
                plot(1:8:n,real(eval(strcat('Phi2',B_cond,'(1:8:n,k)'))), 'ro-', 'markersize',4,'linewidth',1.2);
                set(gca,'Xtick',[0,256,512,768,1024])
            end
            xlim([0,n])
            ylim([-.1,.1])
            h = legend('True', 'DMDc', 'C-DMDc', 'CS-DMDc');
            set(h, 'fontsize', 13, 'location', 'southeast', 'box','on')
            set(gca,'fontsize',15)
            set(gcf,'position',[0, 0,600,150])
            set(gcf,'PaperPositionMode','auto')
            print('-depsc2', '-loose', [figpath, 'phi_',num2str(k), 'B', B_name, '_C', num2str(CType)]);
        end
    end
    
    %% plot B
    if strcmp(B_cond,'2')
        B0 = normalize(B);
        B02 = normalize(Bhat02); [~, Bhat02] = tune_sign(B, Bhat02, 1);
        B12 = normalize(Bhat12); [~, Bbar12] = tune_sign(B, Bhat12, 1);
        B22 = normalize(Bhat22); [~, Bbar22] = tune_sign(B, Bhat22, 1);
        err_Bbar = norm(B-Bhat02);
        err_BY_C = norm(B-Bhat12);
        err_BY_CS = norm(B-Bhat22);
        figure;
        plot(B0,'color',[0.5 0.5 0.5], 'linewidth',8), hold on
        plot(B02,'k-','linewidth',4)
        plot(B12,'b--','linewidth',2);
        plot((1:8:n),B22(1:8:n),'ro-', 'markersize',4,'linewidth',1.2)
        xlim([0,n])
        ylim([-.1,.1])
        h = legend('True', 'DMDc', 'C-DMDc', 'CS-DMDc');
        set(h, 'fontsize', 13, 'location', 'southeast', 'box','on')
        set(gca,'fontsize',15)
        set(gcf,'position',[0, 0,600,150])
        set(gcf,'PaperPositionMode','auto')
        print('-depsc2', '-loose', [figpath, 'B', B_name, '_C', num2str(CType)]);
    end
    
    %% error report
    fileID = fopen([outpath,'error_B', B_name, '_C', num2str(CType), '.txt'], 'w');
    fprintf(fileID,'%s\n','Normalized error of Modes(%):');
    fprintf(fileID,'%20s = %.3d\n','DMDc (B known)', 100*err_Phi_01, 'DMDc (B unknown)', 100*err_Phi_02, ...
        'cDMDc (B known)', 100*err_Phi_11, 'cDMDc (B unknown)', 100*err_Phi_12, ...
        'csDMDc (B known)', 100*err_Phi_21, 'csDMDc (B unknown)', 100*err_Phi_22);
    fprintf(fileID,'%s\n','Normalized error of B(%):');
    fprintf(fileID,'%15s = %.3d\n','DMDc (B unknown)', err_Bbar, ...)
        'cDMDc', err_BY_C, ...
        'csDMDc', err_BY_CS);
    fclose(fileID);
    end
    
end
