function [bigJ,bigX] = construct_block_mat_iv_hi(ndim,varphi,A_tilde,hs,Q,normQ,N)

for j=1:N
    for k=1:N
        [jmat] = calc_j_matrix(ndim,varphi,A_tilde,hs,Q,normQ,j,k);
        [xmat] = calc_x_matrix(ndim,varphi,A_tilde,hs,Q,normQ,j,k);
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        bigJ(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=jmat;
        bigX(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=xmat;
    end
end
end

