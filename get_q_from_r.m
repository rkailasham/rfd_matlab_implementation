function [Qvec] = get_q_from_r(nb,ndim,b2b)
% Returns the connector vector given the matrix
% containing bead-to-bead vectors.

Qvec=zeros((nb-1),ndim);

for j=1:(nb-1)
    Qvec(j,:)=b2b(j,(j+1),:);
end

end

