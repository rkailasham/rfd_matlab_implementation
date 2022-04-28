function [z_mat] = calc_z_mat_proj_op(Q,normQ,j,k)

a=ret_rouse_val(j,k);
dyad=(Q(k,:)'*Q(k,:))./(normQ(k)*normQ(k));
z_mat=a*dyad;
end

