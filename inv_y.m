function [inv_ymat_j] = inv_y(ndim,varphi,Q,normQ,j)
% returns the analytically known inverse of the
% y matrix.
idmat=eye(ndim);
dyad=(Q(j,:)'*Q(j,:))./(normQ(j)*normQ(j));
eps_iv=2.*varphi;
fac_eps=(eps_iv)/(1.+eps_iv);
inv_ymat_j=idmat-((fac_eps)*dyad);
end

