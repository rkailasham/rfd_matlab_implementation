function [deltaR] = modr_b2b(nb,ndim,b2b)
%returns the magnitude of the bead-to-bead vectors
deltaR=zeros(nb,nb);
r12=zeros(1,ndim);
for nu=1:nb
    for mu=1:(nb-1)
        r12(1,:)=b2b(mu,nu,:);
        modr=norm(r12);
        if(modr<=1e-12)
            modr=1e-12;
        end
        deltaR(mu,nu)=modr;
        deltaR(nu,mu)=modr;
    end
end

end

