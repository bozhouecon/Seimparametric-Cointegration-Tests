clc
clear 
rng(1234567)
%%
T      = 250;
reps   = 10000;
C      = [1,0.4,0.8;0.4,1,0.0;0.8,0.0,1];
sigmas = [1 1 1];
Sigma  = sigmas'*sigmas.*C;
f      = 't3'; % 'G' for Gaussian; 't3' for Student-t3; 'skewt4' for Skewed Student-t4
eps    = fn_eps_p3(T,reps,C,sigmas,f);
%%
p      = 3;
r_0    = 1;
r      = 2;
mu     = [2 3 4]';
alpha  = [-0.3 0 0]';
beta   = [1 0 0]';
Pi     = alpha*beta';
Gamma  = [1,1,0;0,1,1;0,0,1];
Gamma1 = Gamma*0.5;
Gamma2 = zeros(3);
Gamma3 = zeros(3);
GAMMA  = [Gamma1,Gamma2,Gamma3];

k_ref  = 1; % lag order specification for reduced-rank regression by 'jcitest'
bw     = eye(p-r_0)*(4/(p-r_0+2)/T)^(2/(p-r_0+4));
% bw     = eye(p-r_0)*(4/(p-r_0+2)/T)^(1/9);
kernel = 'G';  % 'G' for Gaussian kernel; 'logi' for Logistic kernel

Theta  = -0:-10:-10; 
for i = 1:length(Theta)
    i;
    theta   = Theta(i)/T;
    alpha_l = [zeros(p-r,r),theta*ones(p-r,1)]';
    beta_l  = [zeros(p-r,r),ones(p-r,1)]';
    Pi_h    = Pi + alpha_l*beta_l';
    Y       = fn_ECMrnd(mu,Pi_h,GAMMA,eps);
    
    parfor (s = 1:reps,10)
        % s
        %%%%%% Johansen from Matlab
        [h,~,stat,cValue,mles] = jcitest(Y(2:end,:,s),'model','H1*','Lags',k_ref,'display','off');
        alpha_l_hat = null(mles.r1.paramVals.A.');
        Ys(:,:,s) = cumsum(mles.r1.res,1);
        %%%%%% Johansen-phi & Johansen-fhat
        [rej_Johansen_phi(s)] = test_Johansen_phi(Ys(:,:,s)*alpha_l_hat);
        [rej_Johansen_fhat(s)] = test_Johansen_fhat(Ys(:,:,s)*alpha_l_hat,f,bw,kernel);
    end

    power_Johansen_phi(i) = mean(rej_Johansen_phi);
    power_Johansen_fhat(i) = mean(rej_Johansen_fhat);
end

%%
figure;
plot(-Theta,power_Johansen_phi(1:length(Theta)),'--','color','[0 0.45 0.75]');hold on;
plot(-Theta,power_Johansen_fhat(1:length(Theta)),'--','color','[0.85 0.33 0.10]');hold on;
legend('Johansen-$\phi$','Johansen-$\hat{f}$','Location','southeast','interpreter','latex')


