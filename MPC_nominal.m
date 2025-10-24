[A_r_d,B_r_d]=c2d(A_r,B_r,T_s);
x_r=sdpvar(n_r,H_nom+1);
z_r=sdpvar(n_y,H_nom);
u=sdpvar(n_u,H_nom);
cost_z=sdpvar(1,H_nom+1);
cost_u=sdpvar(1,H_nom);
z_ref=z_max*ones(1,H_nom+1);
w_max=0;%no disturbances
%initial state 
obj=R*sum(cost_u)+sum(cost_z);
%constraints=[];
constraints=[cost_u>=0,cost_z>=0];
constraints=[constraints,...
     x_r(:,1)==x_0r];
for k=1:H_nom
%dynamcs
constraints=[constraints,...
    x_r(:,k+1)==A_r_d*x_r(:,k)+B_r_d*u(:,k)];
%input constraints
constraints=[constraints,cost_u(:,k)<=u_max];
%tightened constraints
constraints=[constraints,z_r(:,k)== C_r*x_r(:,k)];
constraints=[constraints,z_r(:,k)<=z_max];
%cost
constraints=[constraints,z_r(:,k)-z_ref(:,k)<=cost_z(k)];
constraints=[constraints,-(z_r(:,k)-z_ref(:,k))<=cost_z(k)];
constraints=[constraints,[eye(n_u),-eye(n_u)]* u(:,k)<=cost_u(k)];
end
opts = sdpsettings('solver', 'mosek');
solve=optimize(constraints,obj,opts);
solve.info
u_ROM=value(u);
z_ROM=value(z_r);
%fprintf('Solve time nominal: %.6f [ms]', solve.solvertime*1e3);
%% plot
%simulate exact
t=linspace(0,T_s*H_nom,H_nom+1);
t=t(1:end-1);
x_full=zeros(n_x,H_nom+1);
[A_d,B_d]=c2d(A,B,T_s);
z_exact=[];
for k=1:H_nom
    x_full(:,k+1)=A_d*x_full(:,k)+B_d*u_ROM(:,k);
    z_exact(:,k)=C*x_full(:,k);
end
% plot
figure(1)
hold on
plot([1,H_nom]*T_s,z_max*[1,1],'black--','linewidth',2)
plot(t,z_ROM,'blue','linewidth',2)
plot(t,z_exact,'red','linewidth',2)   
axis([0,T_s*(H_nom),min([z_exact,z_ROM])*1.1,z_max*1.2])
legend({'Constraint','ROM-prediction','Full-order simulation'},'location','southeast')
xlabel('Time $t$ [s]','interpreter','latex')
ylabel('Position $z [m]$','interpreter','latex')
set(gca, 'fontname','Arial','fontsize',16)
set(gca, 'YTick',[0 0.5 1])
axis([0,20,min([z_exact,z_ROM])*1.1,z_max*1.2])
print('OCP_naive_zoom','-depsc')