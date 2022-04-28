function [div,temp] = div_v_jk_dir_calc_piece_together(i,k,j,m,n,f,s,...
    Q,N,varphi,ndim,del)

%Test function for calculating the divergence of the tensor
%evaluated in "piecing_together.m"



p=(varphi/((2*varphi)+1))^2;
initQ=Q;
temp=zeros(ndim,ndim);


% %Value of "i". The "M" value that is evaluated
% i=2;
% %Value of "k". The "P" value that is evaluated
% k=2;
% %value of "j". Connector vector with respect to which gradient is measured
% j=3;
% 
% %Next we consider the term (D_m/D_n)
% m=5;
% n=3;
% 
% %Next we consider the tensor, T=(\bm{Q}_f\bm{Q}_s)/(Q_fQ_s)
% f=3;
% s=3;

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
        tens_eval_min_in_func=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);
        Q=initQ; %resetting to initial configuration
        %one step forward
        Q(j,col)=Q(j,col)+del;
        normQ = construct_norm(Q,N);
        L = constructL(Q,normQ,N);
        M = fwd_coeff_all(p,L,N);
        P = bkwd_coeff_all(p,L,N);
        tens_eval_plus_in_func=setup_tensor(i,k,j,m,n,f,s,Q,normQ,N,L,p);  
        Q=initQ; %resetting to initial configuration
        temp(row,col)=(tens_eval_plus_in_func(col,row)-tens_eval_min_in_func(col,row));
    end     
end

sumrow=sum(temp,2); %summing along each row
div=(sumrow')./(2*del); %so that it has dimensions of (1,ndim);

% div(1)=temp(1,1)+temp(1,2)+temp(1,3);
% div(2)=temp(2,1)+temp(2,2)+temp(2,3);
% div(3)=temp(3,1)+temp(3,2)+temp(3,3);

end

