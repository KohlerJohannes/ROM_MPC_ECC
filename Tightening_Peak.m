%% peak
r_val=[zeros(n_w,H_validate);xr_val;u_traj_validate];
sys_delta=ss(-lambda_peak,gamma_peak,1,0);
r_val_squared=sum(r_val.^2,1);
delta_peak=lsim(sys_delta,r_val_squared,t_validate,0);
delta_z_peak=sqrt(delta_peak*gamma_peak*lambda_peak);
