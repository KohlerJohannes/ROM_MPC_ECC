%% b) Lohning et al. 2014
%compute P with Lyap using eigenvalue for beta
beta=-max(real(eig(A)))*0.9;%0.9 to leave a small buffer
%
P_lohning=sdpvar(n_x);
con=[A'*P_lohning+P_lohning*A+2*beta*P_lohning+eye(n_x)<=0; P_lohning>=C'*C]
obj=trace(P_lohning);
optimize(con,obj)
%
P_lohning=value(P_lohning);
%scale P s.t. P>=C'*C
alpha=sqrt(max(eig(P_lohning)));%/min(eig(P_lohning)));
%compared to the paper, alpha can possibly reduced if C`\top C<P
r_lohning=[xr_val;u_traj_validate];
r_lohning_norm=sqrt(sum(r_lohning.^2,1));
sys_Lohning=ss(-beta,alpha,1,0);
delta_Lohning=lsim(sys_Lohning,r_lohning_norm,t_validate,0);
 