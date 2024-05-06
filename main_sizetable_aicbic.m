% clc
clear
rng(1234567)
%%
f_set = ["G","t3","skewt4"];
T      = 100;

reps   = 10000;
C      = [1,0.4,0.8;0.4,1,0.0;0.8,0.0,1];
sigmas = [1 1 1];
Sigma  = sigmas'*sigmas.*C;
p      = 3;
r_0    = 1;
r      = 2;
mu     = [2 3 4]';
tau    = [1 1 1]'*0.25;
alpha  = [-0.3 0 0]';
beta   = [1 0 0]';
Pi     = alpha*beta';

Gamma  = [1,1,0;0,1,1;0,0,1];
Gamma1 = Gamma*0.5;
Gamma2 = Gamma*0.3;
Gamma3 = Gamma*0.0;
GAMMA  = [Gamma1,Gamma2,Gamma3];

theta  = 0/T;
alpha_l = [zeros(p-r,r),theta*ones(p-r,1)]';
beta_l  = [zeros(p-r,r),ones(p-r,1)]';
Pi_h    = Pi + alpha_l*beta_l';

bw     = eye(p-r_0)*(4/(p-r_0+2)/T)^(2/(p-r_0+4));
% bw     = eye(p-r_0)*(4/(p-r_0+2)/T)^(1/9);
kernel = 'G';  % G logi

for i = 1:length(f_set)
    f      = f_set(i);
    eps    = fn_eps_p3(T,reps,C,sigmas,f);
    Y      = fn_ECMrnd(mu,Pi_h,GAMMA,eps);
    for k_ref = 1:6
        parfor (s = 1:reps,10)
        % for s = 1:reps
            [h,~,stat,cValue,mles] = jcitest(Y(2:end,:,s),'model','H1*','Lags',k_ref,'display','off');
            [aic(k_ref,s),bic(k_ref,s)] = aicbic(mles.r1.rLL,p*2+p^2*k_ref,T);
            mlesr1(k_ref,s) = mles.r1;
        end
    end
    [~,I_aic] = min(aic,[],1);
    [~,I_bic] = min(bic,[],1);
    parfor (s = 1:reps,10)
    % for s = 1:reps
        alpha_l_hat_aic = null(mlesr1(I_aic(s),s).paramVals.A.');
        alpha_l_hat_bic = null(mlesr1(I_bic(s),s).paramVals.A.');
        Ys_aic = cumsum(mlesr1(I_aic(s),s).res,1);
        Ys_bic = cumsum(mlesr1(I_bic(s),s).res,1);
        %%%%%% Johansen phi & Johansen fhat
        [rej_Joaic_phi(s)] = test_Johansen_phi(Ys_aic*alpha_l_hat_aic);
        [rej_Jobic_phi(s)] = test_Johansen_phi(Ys_bic*alpha_l_hat_bic);
        [rej_Joaic_fhat(s)] = test_Johansen_fhat(Ys_aic*alpha_l_hat_aic,f,bw,kernel);
        [rej_Jobic_fhat(s)] = test_Johansen_fhat(Ys_bic*alpha_l_hat_bic,f,bw,kernel);
    end
    power_Joaic_phi(i) = mean(rej_Joaic_phi);
    power_Jobic_phi(i) = mean(rej_Jobic_phi);
    power_Joaic_fhat(i) = mean(rej_Joaic_fhat);
    power_Jobic_fhat(i) = mean(rej_Jobic_fhat);

    clear Y Ys_aic Ys_bic 

    Y = fn_ECMrnd_trd(mu,tau,Pi_h,GAMMA,eps);
    for k_ref = 1:6
        parfor (s = 1:reps,10)
        % for s = 1:reps
            [h,~,stat,cValue,mles] = jcitest(Y(2:end,:,s),'model','H*','Lags',k_ref,'display','off');
            [aic(k_ref,s),bic(k_ref,s)] = aicbic(mles.r1.rLL,p*2+p^2*k_ref,T);
            mlesr1(k_ref,s) = mles.r1;
        end
    end
    [~,I_aic] = min(aic,[],1);
    [~,I_bic] = min(bic,[],1);
    parfor (s = 1:reps,10)
    % for s = 1:reps
        alpha_l_hat_aic = null(mlesr1(I_aic(s),s).paramVals.A.');
        alpha_l_hat_bic = null(mlesr1(I_bic(s),s).paramVals.A.');
        Ys_aic = cumsum(mlesr1(I_aic(s),s).res,1);
        Ys_bic = cumsum(mlesr1(I_bic(s),s).res,1);
        %%%%%% Johansen phi & Johansen fhat
        [rej_SLaic_phi(s)] = test_SL_phi(Ys_aic*alpha_l_hat_aic);
        [rej_SLbic_phi(s)] = test_SL_phi(Ys_bic*alpha_l_hat_bic);
        [rej_SLaic_fhat(s)] = test_SL_fhat(Ys_aic*alpha_l_hat_aic,f,bw,kernel);
        [rej_SLbic_fhat(s)] = test_SL_fhat(Ys_bic*alpha_l_hat_bic,f,bw,kernel);
    end
    power_SLaic_phi(i) = mean(rej_SLaic_phi);
    power_SLbic_phi(i) = mean(rej_SLbic_phi);
    power_SLaic_fhat(i) = mean(rej_SLaic_fhat);
    power_SLbic_fhat(i) = mean(rej_SLbic_fhat);
    clear Y Ys_aic Ys_bic
end

%%

table_size = table({'$f = Gaussian$';'$f = t_3$';'$f = Skewed-t_4$'},...
    power_Joaic_phi'*100,power_Jobic_phi'*100,power_Joaic_fhat'*100,power_Jobic_fhat'*100, ...
    power_SLaic_phi'*100,power_SLbic_phi'*100,power_SLaic_fhat'*100,power_SLbic_fhat'*100)

% table2latex(table_size,append('table_size_aicbic_T',num2str(T),'.tex'))






