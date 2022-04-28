function [flowMat] = flowMat_eval_direct(varphi,L,Q,normQ,N,ndim)
%Constructing the flow block matrix, \bm{\mathcal{D}}
%using forward and backward substitution coefficients

ndiff=ndim*N;
flowMat=zeros(ndiff,ndiff);
fac=(varphi/((2*varphi)+1));
p=(varphi/((2*varphi)+1))^2;
M = fwd_coeff_all(p,L,N);
P = bkwd_coeff_all(p,L,N);

for j=1:N
    for k=1:N
        [u_hat] = u_hat_jl_calc(j,k,N);        
        [u_numeric] = unwrap_u_direct(u_hat,0,L,varphi,ndim,M,P,Q,normQ);
        idmat=kron_delta(j,k)*eye(ndim); %identity matrix
        a=idmat-(fac*u_numeric);
        drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        flowMat(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=a;
    end
end

end

