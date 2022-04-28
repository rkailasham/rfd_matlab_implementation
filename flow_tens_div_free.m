function [flowMat_div_free] = flow_tens_div_free(N,varphi,Q,normQ,ndim)
%Returns the flow tensor calculated 
%using linear algebra procedure
%[Loads of matrix inversion involved]
%this is specifically for the free-draining case

[bigE,bigV] = construct_block_mat(ndim,varphi,Q,normQ,N);
[bigY] = construct_block_y_mat(ndim,varphi,Q,normQ,N);
% big_Gamma=inv(bigE); %pedagogical step
flowMat_div_free=bigE\bigY; %equivalent to "big_Gamma*bigV" see above
                   %step to understand what is big_Gamma

end

