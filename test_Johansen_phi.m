function [rej] = test_Johansen_phi(Y)

    eps_hat     = diff(Y);
    [T,p]       = size(eps_hat);
    Sigma_hat   = cov(eps_hat);
    We          = [zeros(1,p);cumsum(eps_hat)/sqrt(T)];
    Wv          = We/sqrtm(Sigma_hat);
    CS          = Wv'*[diff(Wv);zeros(1,p)];
    QV          = Wv'*Wv/T;
    LR          = trace(CS'/QV*CS);
    
    rej         = LR>12.3;
end