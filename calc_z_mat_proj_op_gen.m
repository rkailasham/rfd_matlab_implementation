function [z_mat] = calc_z_mat_proj_op_gen(ndim,Q,normQ,A_tilde,j,k)

[A_tilde_jk] = ret_A_tilde_jk(ndim,A_tilde,j,k);
dyad=(Q(k,:)'*Q(k,:))./(normQ(k)*normQ(k));
z_mat=A_tilde_jk*dyad;
end

