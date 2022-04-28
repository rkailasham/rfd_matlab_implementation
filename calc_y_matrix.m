function [ymat_j] = calc_y_matrix(ndim,varphi,Q,normQ,j)
%calculation of Y matrix
idmat=eye(ndim);
[z_mat] = calc_z_mat_proj_op(Q,normQ,j,j);
ymat_j=idmat+(varphi*z_mat);
end

