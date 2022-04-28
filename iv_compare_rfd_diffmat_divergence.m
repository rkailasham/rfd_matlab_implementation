clf;
legend_font_size=20;


%number of springs
N=4;
%number of beads
nb=(N+1);
%number of dimensions
ndim=3;
%size of diffusion matrix
ndsize=N*ndim;
%\varphi=K/zeta
varphi=3;
p=(varphi/((2*varphi)+1));
eps_iv=2*varphi;
%free-draining chains, no HI
hstar=0.0;

div_del=1e-5;

rng('shuffle');


RBead=normrnd(0,1,[nb,ndim]);
[b2b] = b2bvector(nb,ndim,RBead);
Q=get_q_from_r(nb,ndim,b2b);
normQ = construct_norm(Q,N);
L = constructL(Q,normQ,N);

%evaluating the divergence of the diffusion matrix.
%This is a vector with (N) entries, each of which is
%an (ndim)-dimensional vector.

cumul_div=zeros(N,ndim);
div_list=zeros(N,N,ndim);

for k=1:N
    for j=1:N
        [v_hat_jk] = v_hat_jl_calc(j,k,N);
        [div_v_jk_assembly_line] = -p*unwrap_div_v(v_hat_jk,j,0,L,...
                                     varphi,ndim,Q,normQ,N);
        cumul_div(k,:)=cumul_div(k,:)+div_v_jk_assembly_line;  
        div_list(k,j,:)=div_v_jk_assembly_line;
    end
end


%%% Divergence evaluation using the random finite difference approach


rfd_param=1e-5;
% ntraj_list= [10 20 50 100 2e2 5e2 1e3 2e3 5e3 1e4];

ntraj_list= [10 20 50];



[dat,err_dat] = iv_call_rfd_calc_diff_ens_size(varphi,Q,...
                                        N,ndim,rfd_param,ntraj_list);

 %%%%%%%% plotting the results %%%%%%%%%%%%%%
 
 
 %connector vector index to be displayed
k=2;
 
axes1=subplot(3,1,1);
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
axes1.XScale='log';

% xlabel('$N_{\mathrm{trajectories}}$','FontSize',30,'Interpreter','latex');
y1=ylabel('[\boldmath$\nabla\cdot $\textbf{D}$]_{k,x}$','FontSize',30,'Interpreter','latex',...
    'Rotation',90);
% ylim([-0.25 0.1]);
xlim([9. 1e5]);


p2=errorbar(ntraj_list,dat(:,k,1),err_dat(:,k,1),...
    '-d');
p2.MarkerSize=8;
p2.Color=[0.5 0. 0.5];
p2.LineWidth=2;
p2.DisplayName='Random finite difference';
hold on;

h1=refline([0. cumul_div(k,1)]);
h1.LineWidth=3;
h1.Color='k';
h1.DisplayName='Analytical result';
h1.LineStyle='-.';
% h1.HandleVisibility='off';
hold on;
% 
% 
[h,icons,plots,legend_text]=legend({},'Location','northeast','FontSize',legend_font_size,'Interpreter','latex','Box','on');


axes2=subplot(3,1,2);
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
axes2.XScale='log';
% xlabel('$N_{\mathrm{trajectories}}$','FontSize',30,'Interpreter','latex');

y2=ylabel('[\boldmath$\nabla\cdot $\textbf{D}$]_{k,y}$','FontSize',30,'Interpreter','latex',...
    'Rotation',90);% title('$s=3.0$','Interpreter','latex','FontSize',24);
% title('$t/\tilde{\tau}_{1}=3.0$','Interpreter','latex','FontSize',24);

% ylim([-0.1 0.3]);
xlim([9. 1e5]);


p3=errorbar(ntraj_list,dat(:,k,2),err_dat(:,k,2),...
    '-o');
p3.MarkerSize=8;
p3.Color=[0. 0.5 0.];
p3.LineWidth=2;
p3.DisplayName='Random finite difference';
hold on;


h2=refline([0. cumul_div(k,2)]);
h2.LineWidth=1;
h2.Color='r';
h2.DisplayName='Analytical result';
h2.LineStyle='-';
% h2.HandleVisibility='off';


[h,icons,plots,legend_text]=legend({},'Location','southeast','FontSize',legend_font_size,'Interpreter','latex','Box','on');




axes3=subplot(3,1,3);
hold(axes3,'on');
box(axes3,'on');
set(axes3,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
axes3.XScale='log';
xlabel('No. of trajectories','FontSize',30,'Interpreter','latex');
% title('$s=5.0$','Interpreter','latex','FontSize',24);
% title('$t/\tilde{\tau}_{1}=5.0$','Interpreter','latex','FontSize',24);
y3=ylabel('[\boldmath$\nabla\cdot $\textbf{D}$]_{k,z}$','FontSize',30,'Interpreter','latex',...
    'Rotation',90);
% ylim([-0.5 0.2]);
xlim([9. 1e5]);



p4=errorbar(ntraj_list,dat(:,k,3),err_dat(:,k,3),...
    '-s');
p4.MarkerSize=8;
p4.Color=[0.8 0.1 0.8];
p4.LineWidth=2;
p4.DisplayName='Random finite difference';
hold on;

h3=refline([0. cumul_div(k,3)]);
h3.LineWidth=3;
h3.Color='b';
h3.DisplayName='Analytical result';
h3.LineStyle=':';
% h3.HandleVisibility='off';

 % % %common title over all panels
sgtitle(['$N_{\mathrm{b}} =$' num2str(nb) ',\,$k =$' num2str(k),...
    ',$\,\epsilon =$' num2str(eps_iv,'%3.1f') ',\,$h^{*} =$' num2str(hstar,'%3.1f')],'Interpreter',...
    'latex','FontSize',32);


[h,icons,plots,legend_text]=legend({},'Location','best','FontSize',legend_font_size,'Interpreter','latex','Box','on');


 
