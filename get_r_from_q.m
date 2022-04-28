function [RBead] = get_r_from_q(N,ndim,Q)
% Constructs the bead position vectors based on the
% connector vectors. N=number of springs
nb=(N+1);
RBead=zeros(nb,ndim);
for i=1:N
    RBead((i+1),:)=RBead(i,:)+Q(i,:);
end

end

