function epsilon = fn_eps_p3(T,reps,C,sigmas,f)

p = 3;
epsilon = zeros(T,p,reps);

for i = 1:reps  
    if strcmp(f,'G')   
        innovations      = mvnrnd(zeros(p,1),C,T);
    elseif strcmp(f,'t3')
        innovations      = mvtrnd(C,3,T)/sqrt(3/(3-2));
    elseif strcmp(f,'t4')
        innovations      = mvtrnd(C,4,T)/sqrt(4/(4-2));
    elseif strcmp(f,'skewt4')
        U                = copularnd('Gaussian',C,T);
        lambda           = -0.5;
        innovations(:,1) = skewtdis_inv(U(:,1),4,lambda);
        innovations(:,2) = skewtdis_inv(U(:,2),4,lambda);
        innovations(:,3) = skewtdis_inv(U(:,3),4,lambda);
    else
         error('Distribution not yet implemented.')
    end  
    epsilon(:,:,i) = innovations.*repmat(sigmas,T,1);
end
    
end






