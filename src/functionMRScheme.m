function [SE_MR_cent] = functionMRScheme(Hhat,H,D,B,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,R,pilot_stream)

prelogFactor = (1-tau_p/tau_c);

SE_MR_cent = zeros(K,1);


gki_MR = zeros(K,L,K);
gki2_MR = zeros(K,L,K);
Fk_MR = zeros(L,K);

for l = 1:L
    servedUEs = find(D(l,:)==1);
    for ind = 1:length(servedUEs)
        k = servedUEs(ind);
        Fk_MR(l,k) = trace(B(:,:,l,k));
        for i = 1:K
            gki2_MR(i,l,k) = real(trace(B(:,:,l,k)*R(:,:,l,i))); 
            if pilot_stream(k) == pilot_stream(i)
                gki_MR(i,l,k) = real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i)));
                gki2_MR(i,l,k) = gki2_MR(i,l,k) + (real(trace((B(:,:,l,k)/R(:,:,l,k))*R(:,:,l,i))))^2;
            end
        end
    end
end

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
        
        v = Hhatallj_active(:,k);
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        SE_MR_cent(k) = SE_MR_cent(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
    end
end