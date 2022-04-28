function [A_tilde_jk] = ret_A_tilde_jk(ndim,A_tilde,j,k)
%Given the block HI tensor A_tilde of dimensions 3N * 3N
%where N is the number of springs in the chain, function returns
% the 3 x 3 tensor corresponding to the index (j,k)
drow_start_index=((j-1)*ndim)+1;
dcol_start_index=((k-1)*ndim)+1;
drow_fin_index=drow_start_index+(ndim-1);
dcol_fin_index=dcol_start_index+(ndim-1);
A_tilde_jk=A_tilde(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index);
end

