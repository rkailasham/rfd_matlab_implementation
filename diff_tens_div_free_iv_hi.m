function [diffMat_div_free] = diff_tens_div_free_iv_hi(N,varphi,A_tilde,hs,Q,normQ,ndim)
%Returns the diffusion tensor calculated 
%using linear algebra procedure
%[Loads of matrix inversion involved]

[bigJ,bigX] = construct_block_mat_iv_hi(ndim,varphi,A_tilde,hs,Q,normQ,N);
% big_Gamma=inv(bigJ); %pedagogical step
diffMat_div_free=bigJ\bigX; %equivalent to "big_Gamma*bigX" see above
                   %step to understand what is big_Gamma

end

