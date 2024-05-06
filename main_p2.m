% clc
clear 
rng(1234567)
%%
T      = 250;
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
mu     = [0 0]';
Gamma1 = zeros(p);
Gamma2 = zeros(p);
Gamma3 = zeros(p);
GAMMA  = [Gamma1,Gamma2,Gamma3];

aT     = eye(p-r_0)*(4/(p-r_0+2)/T)^(2/(p-r_0+4)); 
% aT     = eye(p-r_0)*(4/(p-r_0+2)/T)^(1/9);
kernel = 'G';  % 'G' for Gaussian kernel; 'logi' for Logistic kernel

Theta  = -0:-10:-10; 

for i = 1:length(Theta)
    i;
    theta   = Theta(i)/T;
    alpha_l = [diag(theta*ones(r,1));zeros(p-r,r)];
    beta_l  = [diag(ones(r,1));ones(p-r,r)];
    Pi_h    = alpha_l*beta_l';
    Y       = fn_ECMrnd(mu,Pi_h,GAMMA,eps);

    parfor (s=1:reps,10)
        %%%%%% Johansen-phi & Johansen-fhat
        [rej_Johansen_phi(s)]        = test_Johansen_phi(Y(:,:,s));
        [rej_Johansen_fhat(s)]       = test_Johansen_fhat(Y(:,:,s),f,aT,kernel);
    end
    
    power_Johansen_phi(i) = mean(rej_Johansen_phi);
    power_Johansen_fhat(i) = mean(rej_Johansen_fhat);
end

%%
figure;
plot(-Theta,power_Johansen_phi(1:length(Theta)),'--','color','[0 0.45 0.75]');hold on;
plot(-Theta,power_Johansen_fhat(1:length(Theta)),'--','color','[0.85 0.33 0.10]');hold on;
legend('Johansen-$\phi$','Johansen-$\hat{f}$','Location','southeast','interpreter','latex')


