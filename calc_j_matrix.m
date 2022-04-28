function [jmat] = calc_j_matrix(ndim,varphi,A_tilde,hs,Q,normQ,j,k)
%calculation of J matrix
idmat=eye(ndim);

if(j==k)
    jmat=idmat;
else
    [z_mat] = calc_z_mat_proj_op_gen(ndim,Q,normQ,A_tilde,j,k);
    ymat_j = ret_ymat(ndim,varphi,hs,Q,normQ,j);
    jmat=varphi*(ymat_j*z_mat);
end

end

