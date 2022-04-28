function [bet] = ret_bet(Qvec,hstar)
% Returns the scalar "beta" which is a function of the
% connector vector
astar=sqrt(pi)*hstar;
normq=norm(Qvec);
alpha=0.75*astar;
[A,B] = ret_rpy_branches(Qvec,hstar);
bet=1.-((alpha./normq)*(A+B));
end

