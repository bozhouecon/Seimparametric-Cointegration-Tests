function [rej] = test_Johansen_fhat(Y,f,aT,kernel)

    eps_hat   = diff(Y);
    [T,p]     = size(eps_hat);
    Sigma_hat = cov(eps_hat);
    v_t       = eps_hat/sqrtm(Sigma_hat);
    v_t_f     = v_t;

    bT = 0.0001;
    
    score_hat = zeros(p,T);
    for t = 1:T
        u = (repmat(v_t(t,:),T,1)-v_t_f)/sqrtm(aT);
        if strcmp(kernel,'G')
            score_hat(:,t) = sum(repmat(prod(exp(-0.5*u.*u),2),1,2).*u/sqrtm(aT),1)/sum(prod(exp(-0.5*u.*u),2)+bT,1);      
        elseif  strcmp(kernel,'logi') 
            score_hat(:,t) = sum(repmat(prod(exp(-u).*(1+exp(-u)).^(-2),2),1,2).*(1+exp(-u)).^(-1).*(1-exp(-u))/sqrtm(aT),1)...
                             /sum(prod(exp(-u).*(1+exp(-u)).^(-2),2)+bT,1);
        end
    end
    
    Wv   = [zeros(1,p);cumsum(v_t)/sqrt(T)];
    Wlfv = cumsum([zeros(1,p);score_hat'])/sqrt(T);
    Blfv = Wlfv - repmat((0:1/T:1)',1,p).*repmat(Wlfv(end,:),T+1,1);
    Jfv  = score_hat*score_hat'/T;
    CS   = Wv'*[diff(Blfv);zeros(1,p)] + mean(Wv)'*Wv(end,:); 
    QV   = (Wv*Jfv)'*Wv/T + (mean(Wv)*(eye(p)-Jfv))'*mean(Wv);
    LR   = trace(CS'/QV*CS);

    
    if strcmp(f,'G')  
        cv = 12.3;     
    elseif strcmp(f,'t3') 
        cv = 11.09; 
    elseif strcmp(f,'skewt4') 
        cv = 11.3; 
    else
        error('Critical value not yet implemented.')
    end
    
    rej            = LR>cv;
end
