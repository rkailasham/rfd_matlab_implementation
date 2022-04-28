function [mu_hat_kl_t] = transpose_mu_hat_kl(mu_hat_kl)

%Re-arranges the information stored in "mu_hat_kl", so that its
%transpose may be assembled.
%This would involve exchanging the 9th and 10th columns of
% "mu_hat_jl", so that the tensor (\bm{Q}_f\bm{Q}_s)/(Q_fQ_s) becomes
%(\bm{Q}_s\bm{Q}_f)/(Q_sQ_f). 

%Once the exchange of the two columns
%has been completed, remember that the prefactor,
%(1/(1-M(k)-P(k)) would depend now on the contents of the 10th column.

%"mu_hat_kl" has 10 columns

mu_hat_kl_t=mu_hat_kl;
[nrow,ncol]=size(mu_hat_kl);
temp1=mu_hat_kl(:,ncol-1);
temp2=mu_hat_kl(:,ncol);

mu_hat_kl_t(:,ncol-1)=temp2;
mu_hat_kl_t(:,ncol)=temp1;


end

