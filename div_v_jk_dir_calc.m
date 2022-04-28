function [div] = div_v_jk_dir_calc(v_hat,j,transp_flag,...
    varphi,ndim,Q,N,del)

%Function calculates the divergence of the "V_jk" (or its transpose,
%as specified) directly.
%"transp_flag" indicates if the input tensor is "v_hat" or
%its transpose.
%If this flag is 0, it means "v_hat" has been input.
%If this flag is 1, it means "v_hat_transpose" has been input.

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
        tens_eval_min=unwrap_v_direct(v_hat,transp_flag,...
            L,varphi,ndim,M,P,Q,normQ);
        Q=initQ; %resetting to initial configuration
        %one step forward
        Q(j,col)=Q(j,col)+del;
        normQ = construct_norm(Q,N);
        L = constructL(Q,normQ,N);
        M = fwd_coeff_all(p,L,N);
        P = bkwd_coeff_all(p,L,N);
        tens_eval_plus=unwrap_v_direct(v_hat,transp_flag,...
            L,varphi,ndim,M,P,Q,normQ);   
        Q=initQ; %resetting to initial configuration
        temp(row,col)=(tens_eval_plus(col,row)-tens_eval_min(col,row))./(2.*del);
    end     
end

sumrow=sum(temp,2); %summing along each row
div=(sumrow'); %so that it has dimensions of (1,ndim);
Q=initQ; %for safety

% div(1)=temp(1,1)+temp(1,2)+temp(1,3);
% div(2)=temp(2,1)+temp(2,2)+temp(2,3);
% div(3)=temp(3,1)+temp(3,2)+temp(3,3);

end

