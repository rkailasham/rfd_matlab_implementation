function [omega_ij] = build_rpy(ndim,hstar,rs,rvec,i,j)

% Returns the dimensionless hydrodynamic 
% interaction tensor, \omega_{ij}
% rvec_{ij} = R_{j}-R_{i} and rs=|rvec|

alpha=0.75*sqrt(pi)*hstar;
omega_ij=zeros(ndim,ndim);

[A,B] = ret_rpy_branches(rvec,hstar);

if(i~=j)
    term1=A*eye(ndim);
    dyad=(rvec'*rvec)./(rs*rs);
    term2=B*dyad;
    omega_ij=(alpha./rs)*(term1+term2);
end


end

