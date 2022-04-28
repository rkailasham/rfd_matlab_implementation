function [xmat] = calc_x_matrix(ndim,varphi,A_tilde,hs,Q,normQ,j,k)
%calculation of X matrix
[A_tilde_jk] = ret_A_tilde_jk(ndim,A_tilde,j,k);
ymat_j = ret_ymat(ndim,varphi,hs,Q,normQ,j);
xmat=ymat_j*A_tilde_jk;
end

