clear
clc

T     = 1000;  
reps  = 50000;
rho   = 0;
p     = 2;
Sigma = [1,0;0,1];
% Jf    = [1,0;0,1];         % Gaussian 
% Jf    = [2.13,0;0,2.13];   % t3 
Jf    = [1.91,0;0,1.91];   % skewt4 
SIGMA = [Sigma,eye(p);eye(p),Jf];

for i = 1:reps
    innovs = mvnrnd(zeros(4,1),SIGMA,T);
    Zey    = [zeros(1,1);cumsum(innovs(:,1))/sqrt(T)]; 
    Zex    = [zeros(1,1);cumsum(innovs(:,2))/sqrt(T)];
    Zlfy   = [zeros(1,1);cumsum(innovs(:,3))/sqrt(T)];   
    Zlfx   = [zeros(1,1);cumsum(innovs(:,4))/sqrt(T)];   
    Ze     = [Zey Zex];
    Zlf    = [Zlfy Zlfx];
    Blf    = Zlf-repmat((0:1/T:1)',1,2).*repmat(Zlf(end,:),T+1,1);
    
    CS     = Ze'*[diff(Zlf);zeros(1,p)]+mean(Ze)'*(Ze(end,:)-Zlf(end,:));  
    QV     = (Ze*Jf)'*Ze/T+(mean(Ze)*(eye(p)-Jf))'*mean(Ze);
    LR(i)  = trace(CS'/QV*CS);
end

LR_cv = quantile(LR,0.95)



    

    
   
  
