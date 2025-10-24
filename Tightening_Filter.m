%this code re-computes tightening based on filter outside of solver as a
%sanity check (should return delta_z_opt)
%1. input
r_val=[zeros(n_w,H_validate);xr_val;u_traj_validate];
%2. filtered input
sys_filt=ss(A_psi,B_psi,C_psi,D_psi);
r_psi_sim=lsim(sys_filt,r_val',t_validate,zeros(n_psi,1))';
r_psi_sim_square=sum(r_psi_sim.^2,1);
sys_delta_filt=ss(-lambda_filtered,gamma_filtered,1,0);
delta_filt=lsim(sys_delta_filt,r_psi_sim_square,t_validate,0);
delta_filt_z=sqrt(delta_filt*gamma_filtered*lambda_filtered);
