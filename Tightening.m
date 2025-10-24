%% Compare size of constraint tightening for different methods
%simulate with fine grid for high accuracy
u_traj_validate=repelem(u_opt,T_s/T_s_validate);
H_validate=length(u_traj_validate);
t_validate=linspace(0,T_s_validate*H_validate,H_validate+1);
t_validate=t_validate(1:end-1);
%simulate ROM
[yr_val,~,xr_val]=lsim(sys_r,u_traj_validate,t_validate,x_0r);
xr_val=xr_val';
%a) Uniform: Lorenzetti et al. 2022 -> z_worstcase
Tightening_uniform;
%b) Lohning et al. 2014-> Delta_lohning (needs to solve LMI)
Tightening_Lohning;
%c) peak -> delta_z_peak
Tightening_Peak
%d) peak with filter -> delta_filt_z (should be same as delta_z_opt)
Tightening_Filter
%% plot tightening
figure(3)
hold on
t_validate(1)=1e-6;z_worstcase(1)=1e-6;delta_Lohning(1)=1e-6;delta_z_peak(1)=1e-6;delta_filt_z(1)=1e-6;%to ensure it can be plotted
plot(t_validate,z_worstcase,'linewidth',2)
plot(t_validate,delta_Lohning,'linewidth',2)
plot(t_validate,delta_z_peak,'linewidth',2)
plot(t_validate,delta_filt_z,'linewidth',2)
plot(t_validate,ones(size(t_validate)),'black--','linewidth',1)
legend({'Uniform','Input-dependent','Peak','Peak-filter','Size of constraints'})%,'location','northwest')
set(gca, 'YScale', 'log')
yticks(10.^[-3 0 3 6])
%set(gca, 'XScale', 'log')
axis([T_s_validate,H*T_s,1e-3,max(delta_Lohning)])
xlabel('Time $t$ [s]','interpreter','latex')
ylabel({'Prediction error $\|e_z\|$ [m]'},'interpreter','latex')
set(gca, 'fontname','Arial','fontsize',16)
print('Tightening_comparison_noLog','-depsc')