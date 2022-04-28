function [grad_val] = derv_inv_mi_pk_calc(i,k,j,p,N,Q,ndim,del)
%This function calculates the gradient of [1/(1-M_i-P_k)] w.r.t \bm{Q}_j
%numerically
initQ=Q;
grad_val=zeros(1,ndim);

for chdim=1:ndim
    %backward step
    Q(j,chdim)=Q(j,chdim)-del;
    fval_min=eval_inv_mi_pk(i,k,p,N,Q);
    Q=initQ; %resetting to initial configuration
    %forward step
    Q(j,chdim)=Q(j,chdim)+del;
    fval_plus=eval_inv_mi_pk(i,k,p,N,Q);
    grad_val(chdim)=(fval_plus-fval_min)./(2*del);
    Q=initQ; %resetting to initial configuration
end

Q=initQ; %for safety


end

