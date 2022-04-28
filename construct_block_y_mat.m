function [bigY] = construct_block_y_mat(ndim,varphi,hs,Q,normQ,N)

for j=1:N
    [ymat_j] = ret_ymat(ndim,varphi,hs,Q,normQ,j);
    for k=1:N
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        if(j==k)
            bigY(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=ymat_j;
        else
            bigY(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=0.;
        end
    end
end
end

