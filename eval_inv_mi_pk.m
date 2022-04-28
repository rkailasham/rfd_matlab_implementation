function [func_val] = eval_inv_mi_pk(i,k,p,N,Q)
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
M_i = fwd_coeff_specific(p,L,i);
P_k = bkwd_coeff_specific(p,L,k,N);
func_val=1./(1.-M_i-P_k);
end

