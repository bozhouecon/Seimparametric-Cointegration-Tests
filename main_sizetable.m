% clc
clear
% rng(1234567)
%%
f_set = ["G","t3","skewt4"];
T      = 250;

reps   = 20000;
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
theta  = 0/T;
Gamma = [1,1,0;0,1,1;0,0,1];
alpha_l = [zeros(p-r,r),theta*ones(p-r,1)]';
beta_l  = [zeros(p-r,r),ones(p-r,1)]';
Pi_h    = Pi + alpha_l*beta_l';

bw     = eye(p-r_0)*(4/(p-r_0+2)/T)^(2/(p-r_0+4));
% bw     = eye(p-r_0)*(4/(p-r_0+2)/T)^(1/9);
kernel = 'G';  % G logi


for i = 1:length(f_set)
    f      = f_set(i);
    eps    = fn_eps_p3(T,reps,C,sigmas,f);

    for k_tru = 1:3
        if k_tru == 1
            Gamma1 = Gamma*0.5;
            Gamma2 = zeros(p);
            Gamma3 = zeros(p);
        elseif k_tru == 2
            Gamma1 = Gamma*0.5;
            Gamma2 = Gamma*0.3;
            Gamma3 = zeros(p);
        elseif k_tru == 3
            Gamma1 = Gamma*0.5;
            Gamma2 = Gamma*0.3;
            Gamma3 = Gamma*0.1;
        end

        GAMMA = [Gamma1,Gamma2,Gamma3];

        for k_ref = 1:3
            Y = fn_ECMrnd(mu,Pi_h,GAMMA,eps);

            parfor (s = 1:reps,10)
                %%%%%% Johansen from Matlab
                [~,~,stat,cValue,mles] = jcitest(Y(2:end,:,s),'model','H1*','Lags',k_ref,'display','off');
                alpha_l_hat = null(mles.r1.paramVals.A.');
                Ys(:,:,s) = cumsum(mles.r1.res,1);
                % aicbic(mles.r1.rLL,p*2+p^2*k_ref,T)
                %%%%%% Johansen phi & Johansen fhat
                % [h,~,~,~] = jcitest(Ys(:,:,s)*alpha_l_hat,'model','H2','display','off');
                % [rej_Jo_phi(s)]      = table2array(h(1,1));
                [rej_Jo_phi(s)] = test_Johansen_phi(Ys(:,:,s)*alpha_l_hat);
                [rej_Jo_fhat(s)] = test_Johansen_fhat(Ys(:,:,s)*alpha_l_hat,f,bw,kernel);
            end
            power_Johansen_phi(k_tru+(i-1)*3,k_ref) = mean(rej_Jo_phi);
            power_Johansen_fhat(k_tru+(i-1)*3,k_ref) = mean(rej_Jo_fhat);
            clear Ys
        end
    end
end
%%

for i = 1:length(f_set)
    f      = f_set(i);
    eps    = fn_eps_p3(T,reps,C,sigmas,f);

    for k_tru = 1:3
        if k_tru == 1
            Gamma1 = Gamma*0.5;
            Gamma2 = zeros(p);
            Gamma3 = zeros(p);
        elseif k_tru == 2
            Gamma1 = Gamma*0.5;
            Gamma2 = Gamma*0.3;
            Gamma3 = zeros(p);
        elseif k_tru == 3
            Gamma1 = Gamma*0.5;
            Gamma2 = Gamma*0.3;
            Gamma3 = Gamma*0.1;
        end

        GAMMA = [Gamma1,Gamma2,Gamma3];

        for k_ref = 1:3

            Y = fn_ECMrnd_trd(mu,tau,Pi_h,GAMMA,eps);

            parfor (s = 1:reps,10)
                %%%%%% Johansen from Matlab
                [~,~,stat,cValue,mles] = jcitest(Y(2:end,:,s),'model','H*','Lags',k_ref,'display','off');
                alpha_l_hat = null(mles.r1.paramVals.A.');
                Ys(:,:,s) = cumsum(mles.r1.res,1);
                % aicbic(mles.r1.rLL,p*2+p^2*k_ref,T)
                %%%%%% Johansen phi & Johansen fhat
                % [h,~,~,~] = jcitest(Ys(:,:,s)*alpha_l_hat,'model','H2','display','off');
                % [rej_Jo_phi(s)]      = table2array(h(1,1));
                [rej_SL_phi(s)] = test_SL_phi(Ys(:,:,s)*alpha_l_hat);
                [rej_SL_fhat(s)] = test_SL_fhat(Ys(:,:,s)*alpha_l_hat,f,bw,kernel);
            end

            power_SL_phi(k_tru+(i-1)*3,k_ref) = mean(rej_SL_phi);
            power_SL_fhat(k_tru+(i-1)*3,k_ref) = mean(rej_SL_fhat);
            clear Ys
        end
    end
end


%%

table_size = table(repmat({'$k=2$';'$k=3$';'$k=4$'},length(f_set),1),...
    power_Johansen_phi(:,1)*100,power_Johansen_phi(:,2)*100,power_Johansen_phi(:,3)*100, ...
    power_Johansen_fhat(:,1)*100,power_Johansen_fhat(:,2)*100,power_Johansen_fhat(:,3)*100, ...
    power_SL_phi(:,1)*100,power_SL_phi(:,2)*100,power_SL_phi(:,3)*100, ...
    power_SL_fhat(:,1)*100,power_SL_fhat(:,2)*100,power_SL_fhat(:,3)*100)

% table2latex(table_size,append('table_size_T',num2str(T),'.tex'))



