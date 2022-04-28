function [div] = div_diffmat_standard_finite_diff(varphi,Q,N,ndim,j,div_del)

% Function calculates the divergence of the expanded diffusion
% matrix using STANDARD FINITE DIFFERENCE.


% "j" refers to the \bm{Q}_j w.r.t which the divergence
% will be taken.
% div_del refers to the dicretization width.

initQ=Q;
temp=zeros(N*ndim,N*ndim);

vecQ=reshape(Q,[N*ndim 1]);
initvecQ=vecQ;


%In MATLAB notation, remember that the first index on A(m,n)
%refers to the row, and the second refers to the column.

%In the following calculation, it is more efficient to traverse row-wise.

for col=1:(N*ndim)
    for row=1:(N*ndim)
        %one step backward
        vecQ(j)=vecQ(j)-div_del;
        normQ = construct_norm(vecQ,N);
        L = constructL(vecQ,normQ,N);
        [Dminus] = diffMat_eval_direct(varphi,L,Q,normQ,N,ndim);
        vecQ=initvecQ; %resetting to initial configuration
        %one step forward
        vecQ(j)=vecQ(j)+div_del;
        normQ = construct_norm(vecQ,N);
        L = constructL(vecQ,normQ,N);
        [Dplus] = diffMat_eval_direct(varphi,L,Q,normQ,N,ndim); 
        vecQ=initvecQ; %resetting to initial configuration
        temp(row,col)=(Dplus(col,row)-Dminus(col,row))./(2.*div_del);
    end     
end

sumrow=sum(temp,2); %summing along each row
div=(sumrow'); %so that it has dimensions of (1,ndim);
Q=initQ; %for safety

% div(1)=temp(1,1)+temp(1,2)+temp(1,3);
% div(2)=temp(2,1)+temp(2,2)+temp(2,3);
% div(3)=temp(3,1)+temp(3,2)+temp(3,3);

end

