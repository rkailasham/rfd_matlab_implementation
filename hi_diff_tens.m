function [a_tilde_expanded] = hi_diff_tens(nb,ndim,hstar,deltaR,b2b)

%Building the expanded version of the block diffusion tensor
%Number of springs, N = (nb-1)
N=(nb-1);
ndiff=N*ndim;

rvec=zeros(1,ndim);

a_tilde_expanded=zeros(ndiff,ndiff);

for k=1:N
    for j=1:N
        a_r_jk=ret_rouse_val(j,k)*eye(ndim);   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        rvec(1,:)=b2b(j,k,:);
        rs=deltaR(j,k);
        temp_jk=build_rpy(ndim,hstar,rs,rvec,j,k);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rvec(1,:)=b2b((j+1),(k+1),:);
        rs=deltaR((j+1),(k+1));              
        temp_j1k1=build_rpy(ndim,hstar,rs,rvec,(j+1),(k+1)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rvec(1,:)=b2b((j),(k+1),:);
        rs=deltaR((j),(k+1)); 
        temp_jk1=build_rpy(ndim,hstar,rs,rvec,(j),(k+1)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rvec(1,:)=b2b((j+1),(k),:);
        rs=deltaR((j+1),(k)); 
        temp_j1k=build_rpy(ndim,hstar,rs,rvec,(j+1),(k));   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_til_jk=(a_r_jk)+temp_jk+temp_j1k1-temp_jk1-temp_j1k;
                drow_start_index=((j-1)*ndim)+1;
        dcol_start_index=((k-1)*ndim)+1;
        drow_fin_index=drow_start_index+(ndim-1);
        dcol_fin_index=dcol_start_index+(ndim-1);
        a_tilde_expanded(drow_start_index:drow_fin_index,dcol_start_index:dcol_fin_index)=a_til_jk;
    end 
end

