function [dat,err_dat] = iv_call_rfd_calc_diff_ens_size(varphi,Q,...
                                        N,ndim,rfd_param,ntraj_list)


%%% Driver routine that calls "rfd_calc.m" and generates
%%% divergences at different ensemble sizes.

size_val=length(ntraj_list);

dat=zeros(size_val,N,ndim);
err_dat=zeros(size_val,N,ndim);


for i=1:length(ntraj_list)
    [divergence,err_divergence] = iv_rfd_calc(varphi,Q,...
        N,ndim,rfd_param,ntraj_list(i));
    dat(i,:,:)=[divergence];
    err_dat(i,:,:)=[err_divergence];
end

end

