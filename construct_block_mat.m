function [bigE,bigV] = construct_block_mat(ndim,varphi,Q,normQ,N)

for j=1:N
    for k=1:N
        [emat] = calc_e_matrix(ndim,varphi,Q,normQ,j,k);
        [vmat] = calc_v_matrix(ndim,varphi,Q,normQ,j,k);
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        bigE(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=emat;
        bigV(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=vmat;
    end
end
end

