function [diffMat_alt] = diffMat_eval_altern(varphi,L,Q,normQ,N,ndim)
%Constructing the diffusion block matrix, \bm{\mathcal{D}}
%but with the order of elements interchanged
ndiff=ndim*N;
diffMat_alt=zeros(ndiff,ndiff);
fac=(varphi/((2*varphi)+1));

for j=1:N
    for k=1:N
        [v_hat] = v_hat_jl_calc(k,j,N);
        [v_numeric] = unwrap_v(v_hat,0,L,varphi,ndim,Q,normQ,N);
        idmat=ret_rouse_val(j,k)*eye(ndim); %identity matrix
        a=idmat-(fac*v_numeric);
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        diffMat_alt(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=a;
    end
end

end

