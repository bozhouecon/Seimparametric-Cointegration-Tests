% clc
clear 
rng(1234567)
%%
T      = 1000;
reps   = 10000;
C      = [1,0;0,1];
sigmas = [1 1];
Sigma  = sigmas'*sigmas.*C;
f      = 'skewt4';
eps    = fn_eps_p2(T,reps,C,sigmas,f);
%%
p      = 2;
r_0    = 0;
r      = 1;
mu     = [3 4]';
tau    = [1 1]'*0.25;
Gamma1 = zeros(p);
Gamma2 = zeros(p);
Gamma3 = zeros(p);
GAMMA  = [Gamma1,Gamma2,Gamma3];

aT     = eye(p-r_0)*(4/(p-r_0+2)/T)^(2/(p-r_0+4));
% aT     = eye(p-r_0)*(4/(p-r_0+2)/T)^(1/9);
kernel = 'G';  % 'G' for Gaussian kernel; 'logi' for Logistic kernel

Theta  = -0:-10:-10;    % -0:-1:-30; 

for i = 1:length(Theta)
    i;
    theta   = Theta(i)/T;
    alpha_l = [zeros(p-r,r),theta*ones(p-r,1)]';
    beta_l  = [zeros(p-r,r),ones(p-r,1)]';
    Pi      = zeros(p);
    Pi_h    = Pi + alpha_l*beta_l';
    Y       = fn_ECMrnd_trd(mu,tau,Pi_h,GAMMA,eps);

    parfor (s=1:reps, 10)
        s;
        %%%%%% SL-phi & SL-fhat
        [rej_SL_phi(s)] = test_SL_phi(Y(:,:,s));
        [rej_SL_fhat(s)] = test_SL_fhat(Y(:,:,s),f,aT,kernel);
    end  
    
    power_SL_phi(i) = mean(rej_SL_phi);
    power_SL_fhat(i) = mean(rej_SL_fhat);
end
% theta_hat = mean(th,1);
%%
figure;
plot(-Theta,power_SL_phi(1:length(Theta)),'--','color','[0 0.45 0.75]');hold on;
plot(-Theta,power_SL_fhat(1:length(Theta)),'--','color','[0.85 0.33 0.10]');hold on;
legend('Johansen-$\phi$','Johansen-$\hat{f}$','Location','southeast','interpreter','latex')





