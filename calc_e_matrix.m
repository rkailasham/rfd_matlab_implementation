function [emat] = calc_e_matrix(ndim,varphi,Q,normQ,j,k)
%calculation of E matrix
idmat=eye(ndim);

if(j==k)
    emat=idmat;
else
    [z_mat] = calc_z_mat_proj_op(Q,normQ,j,k);
    inv_ymat_j = inv_y(ndim,varphi,Q,normQ,j);
    emat=varphi*(inv_ymat_j*z_mat);
end

end

