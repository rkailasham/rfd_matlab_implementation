function [divergence,err_divergence] = iv_rfd_calc(varphi,Q,N,ndim,rfd_param,ntraj)

%%% Implements the Random finite difference technique
%%% for the calulation of the divergence of a 
%%% configuration dependent tensor. See Eq. (26) in
%%% J. Chem. Phys. 140, 134110 (2014) for details.
%%% Link to arxiv version: https://arxiv.org/abs/1401.4198


Q_orig=Q; %saving initial configuration


divergence=zeros(N*ndim,1);
err_divergence=zeros(N*ndim,1);


for i=1:ntraj
    w=normrnd(0,1,[N,ndim]);
    BigW=reshape(w',[],1);
%     BigW=reshape(w,[N*ndim 1]);
    % moving by delta in the +ve direction
    Q=Q+(0.5*rfd_param*w);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    [Dplus] = diffMat_eval_direct(varphi,L,Q,normQ,N,ndim);
    Q=Q_orig;    
    % moving by delta in the -ve direction
    Q=Q-(0.5*rfd_param*w);
    normQ = construct_norm(Q,N);
    L = constructL(Q,normQ,N);
    [Dminus] = diffMat_eval_direct(varphi,L,Q,normQ,N,ndim);
    Q=Q_orig;      
    D_diff=Dplus-Dminus;
    temp_dot=(D_diff*BigW);
    divergence=divergence+(temp_dot);
    err_divergence=err_divergence+(temp_dot.^2);
end

% Methodology described in Oettinger's textbook
% for the calculation of error bars. See solution
% exercise 4.11 on page 311.

divergence=divergence./(ntraj);
err_divergence=err_divergence./(ntraj);
err_divergence=(err_divergence-(divergence.^2))./(ntraj-1);
err_divergence=sqrt(err_divergence);

temp_divergence=divergence./rfd_param;
temp_err_divergence=err_divergence./rfd_param;

mid_div=reshape(temp_divergence,[ndim N]);
mid_err=reshape(temp_err_divergence,[ndim N]);

divergence=mid_div';
err_divergence=mid_err';

% divergence=temp_divergence;
% err_divergence=temp_err_divergence;
% % size(temp_divergence)
% % size(divergence)

end