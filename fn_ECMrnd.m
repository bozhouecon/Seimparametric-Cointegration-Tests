function   X = fn_ECMrnd(mu,Pi,GAMMA,epsilon)
% X = ECMrnd(T,mu,Pi,epsilon)
% This function can be used to simulate a ECM (without lags) starting at 0
% INPUT:
%   * T: number of observations (process is started at 0)
%   * mu: n x 1 vector of intercepts
%   * Pi:   n x n matrix
%   * epsilon:  n x T x reps matrix with innovations
% OUTPUT:
%   The function simulates trajectories X_t according to the error correction model
%   \Delta X_t = mu + Pi X_{t-1} + \epsilon_t
%   dimension of X is n x T x reps 
%% determine dimensions
[T,p,reps] = size(epsilon);
Gamma1 = GAMMA(:,1:p);
Gamma2 = GAMMA(:,(p+1):2*p);
Gamma3 = GAMMA(:,(2*p+1):3*p);
%% simulation process
X = zeros(T+1,p,reps);
for t = 5:T+1
    H  = squeeze(X(t-1,:,:));
    P1 = squeeze(X(t-1,:,:)-X(t-2,:,:));
    P2 = squeeze(X(t-2,:,:)-X(t-3,:,:));
    P3 = squeeze(X(t-3,:,:)-X(t-4,:,:));
    X(t,:,:) = H + Pi*H + Gamma1*P1 + Gamma2*P2 + Gamma3*P3 + squeeze(epsilon(t-1,:,:));
end
X = X(1:T+1,:,:) + repmat(mu',T+1,1,reps); % remove starting value
end