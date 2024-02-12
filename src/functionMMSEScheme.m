function [SE_MMSE] = functionMMSEScheme(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,p)

prelogFactor = (1-tau_p/tau_c);

SE_MMSE =  zeros(K,1);

for n = 1:nbrOfRealizations
    for k = 1:K
        servingAPs = find(D(:,k)==1); %cell-free setup
        La = length(servingAPs);
        
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        Hallj_active = zeros(N*La,K);
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
        v = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        SE_MMSE(k) = SE_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
    end
end