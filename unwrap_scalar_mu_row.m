function [mu_scalar_numeric] = unwrap_scalar_mu_row(mu_hat_rownum,transp_flag,...
    L,p,iv_pf,N)

%Unwraps the scalar part of a given row (rownum) of matrix "mu" and 
%returns its numeric value
%for a given chain configuration. The matrix "mu" is obtained
%as an output of either "mu_hat_kl_calc" or "transpose_mu_hat_kl"
%"transp_flag" indicates if the input tensor is "u_hat" or
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

[nrow,ncol]=size(mu_hat_rownum);

if(~transp_flag)
    kcol=ncol-1;
else
    kcol=ncol;
end

multi_pf=mu_hat_rownum(1)*mu_hat_rownum(2);
ivfac=(iv_pf)^(mu_hat_rownum(3));

if(mu_hat_rownum(4)==1)
    ratio_seq=fwd_coeff_poly(mu_hat_rownum(5),L,p,N)/fwd_coeff_poly(mu_hat_rownum(6),L,p,N);
elseif(mu_hat_rownum(4)==2)
    ratio_seq=bkwd_coeff_poly(mu_hat_rownum(5),L,p,N)/bkwd_coeff_poly(mu_hat_rownum(6),L,p,N);
end
lseq=cont_prod(mu_hat_rownum(7),mu_hat_rownum(8),L);
k=mu_hat_rownum(kcol);
pf=calc_inv_pi_mk(k,k,L,p,N);

mu_scalar_numeric=multi_pf*pf*ivfac*ratio_seq*lseq;

end
