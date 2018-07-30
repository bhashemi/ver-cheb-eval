clc, clear all
% Illustrate linear complexity of three algorithms in the repository.

%% First choose approprate values of lambda so that the corresponding 
% randnfun has appropriate degrees
lam = [0.100000000000001; % corresponds to length 102
0.04; % 225
0.0164  % 453
0.00771  % 901
0.0037             % 1805
0.00255            % 2577
0.00130]            % 4943
sz = numel(lam);
n = zeros(sz,1);
t_LohnerQR = zeros(sz,1);
t_para = t_LohnerQR;
t_ICA_eig_err = t_LohnerQR;

numpts = 20;
rng(1), x = 2*rand(numpts,1)-1; % bary is the best with this choice of points
[~,ind] = sort(x);
x = x(ind);

%%
p_LohnerQR = intval(zeros(numpts,sz));
p_para = p_LohnerQR;
p_ICA_eig_err = p_LohnerQR;

for i=1:sz
    i
    rng(1), f = randnfun(lam(i));
    intCoeffs = f.coeffs;
    n(i) = numel(intCoeffs)
    [p_LohnerQR(:,i), t_LohnerQR(i)] = ICA_QR_err_vec(x, intCoeffs);
    [p_para(:,i), t_para(i)] = ICA_para_err_vec(x, intCoeffs);
    [p_ICA_eig_err(:,i), t_ICA_eig_err(i)] = ICA_eig_err_vec(x, intCoeffs);
end

%%
close all
MS = 'markersize'; ms = 10;
LW = 'linewidth'; lw = 4;
FS = 'fontsize'; fs = 16;

plot(n,n/n(1)*t_LohnerQR(1), 'b:')
hold on
plot(n,t_LohnerQR, 'b-', LW, lw);

plot(n,n/n(1)*t_para(1), 'r-.')
plot(n,t_para, 'r-', LW, lw);

plot(n,n/n(1)*t_ICA_eig_err(1), 'g--')
plot(n,t_ICA_eig_err, 'g-', LW, lw);

leg = legend('t = \alpha n', 'ICA-QR-err', 't = \beta n', 'para', ...
    't = \gamma n', 'ICA-eig-err', 'location', 'nw');
set(leg, FS, fs)
xlabel('degree n', FS, fs+3)
ylabel('time t', FS, fs+3)
title('Verified evaluation at 20 points', FS, fs)
grid on

%%
[mean(-log10(rad(p_para)))'   mean(-log10(rad(p_LohnerQR)))' ...
    mean(-log10(rad(p_ICA_eig_err)))']
