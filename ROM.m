% reduced order
%modal decomposition
[V,Lambda]=eig(A);
%sort eigenvalues
[d,ind] = sort(diag(Lambda),'descend');
V=V(:,ind);
W=inv(V);
V=V(:,end-n_r+1:end);
lambda=d(end-n_r+1:end);
%need to treat real and complex eigenvalues seperately
V_temp=[];
while length(lambda)>0
            if imag(lambda(1))==0
               V_temp=[V_temp,V(:,1)];
               lambda=lambda(2:end);
               V=V(:,2:end);
            elseif imag(lambda(1)+lambda(2))==0
               lambda=lambda(3:end);
                V_temp=[V_temp,real(V(:,1)),imag(V(:,1))];
                V=V(:,3:end);
            else
                error('not implemented')            
            end
end
%add lumped model to ensure offset-free
V=[V_temp,kron(ones(N,1)/sqrt(N),eye(n_i))];
n_r=n_r+2 %number of reduced-order states (6 from slowest eigenvalues and 2 from lumped model)
W=pinv(V)';
     
A_r=W'*A*V;
B_r=W'*B;
C_r=C*V;
sys_r = ss(A_r,B_r,C_r,0);
% error
B_e=[E,(eye(n_x)-V*W')*[A*V B]];
n_input_error=n_w+n_r+n_u;
sys_error=ss(A,B_e,C,zeros(1,n_input_error));