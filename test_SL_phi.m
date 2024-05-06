function [rej] = test_SL_phi(Y)

    eps_hat     = diff(Y);
    [T,p]       = size(eps_hat);
    Sigma_hat   = cov(eps_hat);
    We          = [zeros(1,p);cumsum(eps_hat)/sqrt(T)];
    Wv          = We/sqrtm(Sigma_hat);
    Bv          = Wv - repmat((0:1/T:1)',1,p).*repmat(Wv(end,:),T+1,1);
    CS          = Bv'*[diff(Bv);zeros(1,p)];
    QV          = Bv'*Bv/T;
    LR          = trace(CS'/QV*CS);
    
    rej         = LR>15.92;
end