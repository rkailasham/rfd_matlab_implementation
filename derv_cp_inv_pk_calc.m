function [grad_val] = derv_cp_inv_pk_calc(start,finish,j,p,N,Q,ndim,del)
%This function calculates the gradient of 
%\prod_{i=start}^{i=finish}[1/(1-P_k)] w.r.t \bm{Q}_j

initQ=Q;
grad_val=zeros(1,ndim);

if(j<start)
    grad_val=grad_val;
else
    for chdim=1:ndim
        %backward step
        Q(j,chdim)=Q(j,chdim)-del;
        fval_min=eval_cp_inv_pk(start,finish,p,N,Q);
        Q=initQ; %resetting to initial configuration
        %forward step
        Q(j,chdim)=Q(j,chdim)+del;
        fval_plus=eval_cp_inv_pk(start,finish,p,N,Q);
        grad_val(chdim)=(fval_plus-fval_min)./(2*del);
        Q=initQ; %resetting to initial configuration
    end

end

Q=initQ; %for safety


end

