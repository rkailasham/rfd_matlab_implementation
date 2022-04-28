function [vmat] = calc_v_matrix(ndim,varphi,Q,normQ,j,k)
%calculation of V matrix
idmat=eye(ndim);
a=ret_rouse_val(j,k)*idmat;
inv_ymat_j = inv_y(ndim,varphi,Q,normQ,j);
vmat=inv_ymat_j*a;
end

