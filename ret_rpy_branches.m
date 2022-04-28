function [A,B] = ret_rpy_branches(Qvec,hstar)
% Returns the branches of the RPY tensor for a given value
% of the connector vector and the hydrodynamic interaction 
% parameter, hstar.
astar=sqrt(pi)*hstar;
normq=norm(Qvec);

if(normq>=(2*astar))
    A=1.+((2./3)*(astar/normq)*(astar/normq));
    B=1.-(2.*(astar/normq)*(astar/normq));
else
    A=((4./3)*(normq/astar))-((3./8)*(normq/astar)*(normq/astar));
    B=(1./8)*(normq/astar)*(normq/astar);
end

end

