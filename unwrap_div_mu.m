function [div_mu_numeric] = unwrap_div_mu(mu_hat,j,transp_flag,L,varphi,ndim,Q,normQ,N)

%Takes "mu_hat" or its transpose as input, unwraps it,
%and returns its divergence w.r.t \bm{Q}_j as the numeric output for 
%for a given chain configuration. 
%"transp_flag" indicates if the input tensor is "mu_hat" or
%its transpose.
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

p=(varphi/((2*varphi)+1))^2;
iv_pf=(varphi/((2*varphi)+1));
div_mu_numeric=zeros(1,ndim);

[nrow,ncol]=size(mu_hat);

for r=1:nrow
    [div_mu_row_numeric] = unwrap_div_mu_rowwise(mu_hat(r,:),transp_flag,...
    j,L,Q,normQ,p,iv_pf,N,ndim);
    div_mu_numeric=div_mu_numeric+div_mu_row_numeric;
end


end

