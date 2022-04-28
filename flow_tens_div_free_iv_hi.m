function [flowMat_div_free] = flow_tens_div_free_iv_hi(N,varphi,A_tilde,hs,Q,normQ,ndim)
%Returns the flow tensor calculated 
%using linear algebra procedure
%[Loads of matrix inversion involved]

[bigJ,bigX] = construct_block_mat_iv_hi(ndim,varphi,A_tilde,hs,Q,normQ,N);
[bigY] = construct_block_y_mat(ndim,varphi,hs,Q,normQ,N);
% big_Gamma=inv(bigE); %pedagogical step
flowMat_div_free=bigJ\bigY; %equivalent to "big_Gamma*bigV" see above
                   %step to understand what is big_Gamma

end

