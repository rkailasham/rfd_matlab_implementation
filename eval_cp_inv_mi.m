function [func_val] = eval_cp_inv_mi(start,finish,p,N,Q)
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
M = fwd_coeff_all(p,L,N);
func_val=cont_prod_inv_series(start,finish,M);
end

