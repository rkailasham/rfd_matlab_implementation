function [div,temp] = div_v_jk_dir_calc_dyad(f,s,j,ndim,Q,N,del)

%Function calculates the divergence of the tensor given by 
%T=(\bm{Q}_f\bm{Q}_s)/(Q_fQ_s)
%w.r.t \bm{Q}_j

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
        tens_eval_min=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));
        Q=initQ; %resetting to initial configuration
        %one step forward
        Q(j,col)=Q(j,col)+del;
        normQ = construct_norm(Q,N);
        tens_eval_plus=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));  
        Q=initQ; %resetting to initial configuration
        temp(row,col)=(tens_eval_plus(col,row)-tens_eval_min(col,row))./(2.*del);
    end     
end

sumrow=sum(temp,2); %summing along each row
div=(sumrow'); %so that it has dimensions of (1,ndim);

% div(1)=temp(1,1)+temp(1,2)+temp(1,3);
% div(2)=temp(2,1)+temp(2,2)+temp(2,3);
% div(3)=temp(3,1)+temp(3,2)+temp(3,3);

end

