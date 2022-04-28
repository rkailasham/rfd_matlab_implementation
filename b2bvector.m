function [b2b] = b2bvector(nb,ndim,RBead)
%Computes the bead-to-bead vector between all bead pairs
%The convention is r_{\mu nu}=r_{\nu}-r_{\mu}
b2b=zeros(nb,nb,ndim);
for nu=2:nb
    for mu=1:(nb-1)
        b2b(mu,nu,:) = RBead(nu,:) - RBead(mu,:);
        b2b(nu,mu,:) = -b2b(mu,nu,:);
    end
end

end

