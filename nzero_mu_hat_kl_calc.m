function [func_val] = nzero_mu_hat_kl_calc(mu_hat_kl,N)
%Returns only those rows of "mu_hat_kl"
%which have a non-zero ret_value.

%remember that the input size of mu_hat_kl is (4,10)
%Nonetheless, we do
[nrow,ncol]=size(mu_hat_kl);
zero_ret_val_count=0;

for row=1:nrow
    chk=qlim_flag_chk(mu_hat_kl(row,(ncol-1)),N);
    if(~chk||(mu_hat_kl(row,2)==0))
        mu_hat_kl(row,2)=0; %set ret_val to 0 if pre-multiplying
                           %vector is not legal.
        zero_ret_val_count=zero_ret_val_count+1;                   
    end
end     

tempmat=mu_hat_kl;
func_val=zeros((nrow-zero_ret_val_count),ncol);
build_count=0;

for row=1:nrow
    if(tempmat(row,2))
        build_count=build_count+1;
        func_val(build_count,:)=tempmat(row,:);
    end
end

end

