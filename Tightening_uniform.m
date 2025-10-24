%% Compute robust prediction based on []
%1. determine maximal value of x_r,u
%start at 0 for simplicity
[A_rd_val,B_rd_val]=c2d(A_r,B_r,T_s_validate);
I=eye(n_r);
% Preallocate
temp = zeros(n_r, H_validate);
% Compute all A^k * B in advance
H_max=3*10^3/T_s_validate;%this needs to be long enough, see error check below
AB = zeros(size(A_rd_val,1), H_max); % size [n_x, H_max]
AB(:,1) = B_rd_val;
for k = 2:H_max
    AB(:,k) = A_rd_val * AB(:,k-1);
end
% Multiply by I
temp = I * AB;   % [n_r x H_max]
% Max state
x_max = sum(abs(temp), 2) * u_max;
%this only works for scalar input; otherwise, each entry requires solving a LP
if norm(temp(:,end))>1e-6
  warning('Uniform bound not correctly computed')
  norm(temp(:,end))
  H_max
  %we compute x_max by just taking the maximal possible state for a long horizon;
  %Length of this horizon needs to be chosen problem specific  
end
clear H_max
%
%2. determine maximal value of r
for i=1:n_x
[a,b]=linprog(B_e(i,:),[],[],[],[],-[x_max;u_max],[x_max;u_max]);
d_max(i)=-b;
end
d_max=round(d_max(:));%otherwise, problematic with zero entries
%3. determine maximal error z_e=C*e based on d_max
A_d=expm(A*T_s_validate);

temp = zeros(1,H_validate);
A_powers = cell(1,H_validate);
A_powers{1} = eye(size(A_d));
for k = 2:H_validate
    A_powers{k} = A_d * A_powers{k-1};
end

for k = 1:H_validate
    [~,b] = linprog(C*A_powers{k}, [], [], [], [], -d_max, d_max);
    temp(H_validate-k+1) = -b;
end
z_worstcase=cumsum(flip(abs(temp)));
%sum up to get maximal; flipping needed to get it for every future prediction time