function [divergence,err_divergence] = hi_rfd_calc(nb,ndim,hstar,Q,rfd_param,ntraj)

%%% Implements the Random finite difference technique
%%% for the calulation of the divergence of a 
%%% configuration dependent tensor. See Eq. (26) in
%%% J. Chem. Phys. 140, 134110 (2014) for details.
%%% Link to arxiv version: https://arxiv.org/abs/1401.4198

N=(nb-1);
Q_orig=Q; %saving initial configuration


divergence=zeros(N*ndim,1);
err_divergence=zeros(N*ndim,1);


for i=1:ntraj
    w=normrnd(0,1,[N,ndim]);
    BigW=reshape(w',[],1);
    % moving by delta in the +ve direction
    Q=Q+(0.5*rfd_param*w);
    [RBead] = get_r_from_q(N,ndim,Q);
    [b2b] = b2bvector(nb,ndim,RBead);
    [deltaR] = modr_b2b(nb,ndim,b2b);
    [Dplus] = hi_diff_tens(nb,ndim,hstar,deltaR,b2b);
    Q=Q_orig;    
    % moving by delta in the -ve direction
    Q=Q-(0.5*rfd_param*w);
    [RBead] = get_r_from_q(N,ndim,Q);
    [b2b] = b2bvector(nb,ndim,RBead);
    [deltaR] = modr_b2b(nb,ndim,b2b);
    [Dminus] = hi_diff_tens(nb,ndim,hstar,deltaR,b2b);
    %resetting for safety 
    Q=Q_orig;    
    RBead=0.;
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

end