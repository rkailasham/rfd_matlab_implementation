function [mu_numeric] = unwrap_mu_direct(mu_hat,transp_flag,L,varphi,ndim,M,P,Q,normQ)

%Unwraps the matrix "mu" and returns its numeric value
%for a given chain configuration. The matrix "mu" is obtained
%as an output of either "mu_hat_kl_calc" or "transpose_mu_hat_kl"
%"transp_flag" indicates if the input tensor is "mu_hat" or
%its transpose.
%The suffix "direct" indicates that all quantities are evaluated
%in terms of "M" and "P", without recourse to the polynomial 
%expressions.
%If this flag is 0, it means "mu_hat" has been input.
%If this flag is 1, it means "mu_hat_transpose" has been input.

%%% Here is the structure of "mu" which has 10 columns
%%% Column 1  : Prefactor
%%% Column 2  : ret_val from legality of \bm{Q}_l
%%% Column 3  : iv_exp, exponent to iv prefactor
%%% Column 4  : seq_type (1) for "I", (2) for "D"
%%% Column 5  : Numerator subscript
%%% Column 6  : Denominator subscript
%%% Column 7  : lsq_start, starting index for continued product of L_i's
%%% Column 8  : lsq_fin, finishing index for continued product of L_i's
%%% Column 9  : first index, "f", and 
%%% Column 10 : second index,"s"
%%%             of the dyad (\bm{Q}_f\bm{Q}_s)/(Q_fQ_s)


iv_pf=(varphi/((2*varphi)+1));
mu_numeric=zeros(ndim,ndim);

[nrow,ncol]=size(mu_hat);

for r=1:nrow
    mu_scalar_numeric=unwrap_scalar_mu_row_direct(mu_hat(r,:),...
        transp_flag,L,M,P,iv_pf);
    f=mu_hat(r,(ncol-1));
    s=mu_hat(r,(ncol));
    dyad=(Q(f,:)'*Q(s,:))./(normQ(f)*normQ(s));
    term=mu_scalar_numeric*dyad;
    mu_numeric=mu_numeric+term;
end

end

