function [func_val] = eval_cp_inv_pk(start,finish,p,N,Q)
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);
P = bkwd_coeff_all(p,L,N);
func_val=cont_prod_inv_series(start,finish,P);
end

