function [div_mu_row_numeric] = unwrap_div_mu_rowwise(mu_hat_rownum,transp_flag,...
    j,L,Q,normQ,p,iv_pf,N,ndim)

%Unwraps the divergence of a given row (rownum) of matrix "mu" and 
%returns its numeric value
%for a given chain configuration. 

%The formula is div(fT)=(grad(f)).T + f(div(T)) where "f"
%is scalar and "T" is a tensor

%This is evaluated by calculating the two terms on the RHS of
%the above equation separately, with "fac1=(grad(f)).T" and 
%"fac2=f(div(T))".


%The matrix "mu" is obtained
%as an output of either "mu_hat_kl_calc" or "transpose_mu_hat_kl"
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


div_mu_row_numeric=zeros(1,ndim);

[nrow,ncol]=size(mu_hat_rownum);

if(~transp_flag)
    kcol=ncol-1;
else
    kcol=ncol;
end

%Calculation of fac1=(grad(f)).T

multi_pf=mu_hat_rownum(1)*mu_hat_rownum(2);
ivfac=(iv_pf)^(mu_hat_rownum(3));
pf1=multi_pf*ivfac;

k=mu_hat_rownum(kcol);
[derv1,denom] = derv_pi_mk(k,k,j,p,L,Q,normQ,N);
piece1=(1./denom);

if(mu_hat_rownum(4)==1)
    piece2=fwd_coeff_poly(mu_hat_rownum(5),L,p,N)/fwd_coeff_poly(mu_hat_rownum(6),L,p,N);
    derv2=derv_ratio_i_m_n_q_j(mu_hat_rownum(5),mu_hat_rownum(6),j,N,L,p,Q,normQ);
elseif(mu_hat_rownum(4)==2)
    piece2=bkwd_coeff_poly(mu_hat_rownum(5),L,p,N)/bkwd_coeff_poly(mu_hat_rownum(6),L,p,N);
    derv2=derv_ratio_d_m_n_q_j(mu_hat_rownum(5),mu_hat_rownum(6),j,N,L,p,Q,normQ);
end

[pf2,piece3,derv3] = lseq_partition(j,mu_hat_rownum(7),mu_hat_rownum(8),L,Q,normQ,N);
prefac=pf1*pf2;
first=mu_hat_rownum(9);
sec=mu_hat_rownum(10);
dyad=(Q(first,:)'*Q(sec,:))./(normQ(first)*normQ(sec));

term1=piece1*piece2*derv3;
term2=piece1*piece3*derv2;
term3=piece2*piece3*derv1;

%(grad(f)).T

fac1=prefac*(term1+term2+term3)*dyad;

%Calculation of fac2=f(div(T))

%under construction
f=unwrap_scalar_mu_row(mu_hat_rownum,transp_flag,...
    L,p,iv_pf,N);
divT=derv_q_tens(first,sec,j,Q,normQ);

fac2=f*divT;

div_mu_row_numeric=fac1+fac2;

end



