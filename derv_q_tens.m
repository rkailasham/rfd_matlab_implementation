function [divQ] = derv_q_tens(i,k,j,Q,normQ)
%returns the divergence of [\bm{Q}_i\bm{Q}_k/(Q_iQ_k)]
%w.r.t \bm{Q}_j

if(j==i&&j~=k)
    divQ=(2./(normQ(k)*normQ(j)))*Q(k,:);
elseif(j==k&&j~=i)
    den=normQ(j)*normQ(j)*normQ(j);
    term1=(1./(normQ(i)*normQ(j)))*Q(i,:);
    fac1=(dot(Q(i,:),Q(j,:)))/(normQ(i)*den);
    term2=fac1*Q(j,:);    
    divQ=term1-term2;
elseif(j==i&&j==k)
    divQ=(2./(normQ(j)*normQ(j)))*Q(j,:);
else
    divQ=0;
end

end

