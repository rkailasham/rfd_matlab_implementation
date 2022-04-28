function [ymat_j] = ret_ymat(ndim,varphi,hs,Q,normQ,j)
idmat=eye(ndim);
dyad=(Q(j,:)'*Q(j,:))./(normQ(j)*normQ(j));
eps_iv=2.*varphi;
beta=ret_bet(Q(j,:),hs);
fac_eps=(eps_iv*beta)/(1.+(eps_iv*beta));
ymat_j=idmat-((fac_eps)*dyad);
end

