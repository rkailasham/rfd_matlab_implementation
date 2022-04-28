function [div] = div_mu_kl_dir_calc(mu_hat,j,transp_flag,...
    varphi,ndim,Q,N,del)

%Function calculates the divergence of the "mu_kl" (or its transpose,
%as specified) directly.
%"transp_flag" indicates if the input tensor is "mu_hat" or
%its transpose.
%If this flag is 0, it means "mu_hat" has been input.
%If this flag is 1, it means "mu_hat_transpose" has been input.

%"j" refers to the \bm{Q}_j w.r.t which the divergence
%will be taken.
%del refers to the dicretization width.

p=(varphi/((2*varphi)+1))^2;
initQ=Q;
temp=zeros(ndim,ndim);



%In MATLAB notation, remember that the first index on A(m,n)
%refers to the row, and the second refers to the column.

%In the following calculation, it is more efficient to traverse row-wise.

for col=1:ndim
    for row=1:ndim
        %one step backward
        Q(j,col)=Q(j,col)-del;
        normQ = construct_norm(Q,N);
        L = constructL(Q,normQ,N);
        M = fwd_coeff_all(p,L,N);
        P = bkwd_coeff_all(p,L,N);
        tens_eval_min=unwrap_mu_direct(mu_hat,transp_flag,...
            L,varphi,ndim,M,P,Q,normQ);
        Q=initQ; %resetting to initial configuration
        %one step forward
        Q(j,col)=Q(j,col)+del;
        normQ = construct_norm(Q,N);
        L = constructL(Q,normQ,N);
        M = fwd_coeff_all(p,L,N);
        P = bkwd_coeff_all(p,L,N);
        tens_eval_plus=unwrap_mu_direct(mu_hat,transp_flag,...
            L,varphi,ndim,M,P,Q,normQ);   
        Q=initQ; %resetting to initial configuration
        temp(row,col)=(tens_eval_plus(col,row)-tens_eval_min(col,row))./(2.*del);
    end     
end

sumrow=sum(temp,2); %summing along each row
div=(sumrow'); %so that it has dimensions of (1,ndim);
Q=initQ; %for safety

end
