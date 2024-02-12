function [H_temp,H,B,C] = functionChannelEstimates(R,iterations,L,K,N,tau_p,pilot_stream,p)

H = (randn(L*N,iterations,K)+1i*randn(L*N,iterations,K));

for l = 1:L
    for k = 1:K
        Rsqrt = sqrtm(R(:,:,l,k));
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*H((l-1)*N+1:l*N,:,k);
    end
end

eyeN = eye(N);

Np = sqrt(0.5)*(randn(N,iterations,L,tau_p) + 1i*randn(N,iterations,L,tau_p));
H_temp = zeros(L*N,iterations,K);

if nargout>2
    B = zeros(size(R));
end

if nargout>3
    C = zeros(size(R));
end

for l = 1:L
    for t = 1:tau_p
        yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilot_stream),3) + sqrt(tau_p)*Np(:,:,l,t);
        PsiInv = (p*tau_p*sum(R(:,:,l,t==pilot_stream),4) + eyeN);
        
        for k = find(t==pilot_stream)'
            RPsi = R(:,:,l,k) / PsiInv;
            H_temp((l-1)*N+1:l*N,:,k) = sqrt(p)*RPsi*yp;
            
            if nargout>2
                B(:,:,l,k) = p*tau_p*RPsi*R(:,:,l,k);
            end
            if nargout>3
                C(:,:,l,k) = R(:,:,l,k) - B(:,:,l,k);
            end
        end
    end
end
