function [diffMat_div_free] = diff_tens_div_free(N,varphi,Q,normQ,ndim)
%Returns the diffusion tensor calculated 
%using linear algebra procedure
%[Loads of matrix inversion involved]
%this is specifically for the free-draining case

[bigE,bigV] = construct_block_mat(ndim,varphi,Q,normQ,N);
% big_Gamma=inv(bigE); %pedagogical step
diffMat_div_free=bigE\bigV; %equivalent to "big_Gamma*bigV" see above
                   %step to understand what is big_Gamma

end

