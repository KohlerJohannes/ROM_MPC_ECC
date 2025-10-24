%define dynamics
n_r=6;
n_i=2;
N=50;
if N<3
   error('not programmed') 
end
n_x=n_i*N;
m=1;k=10;d=2*10;
%critical damping: m=1: d=k^2/4
A=zeros(n_x);
%i=1:also contected to ground
A(1,:)=[0,1,zeros(1,n_x-n_i)];
A(2,:)=1/m*[-k-k,-d-d,k,d,zeros(1,n_x-2*n_i)];
for i=2:N-1
A((i-1)*n_i+1,:)=[zeros(1,(i-1)*n_i+1),1,zeros(1,n_x-i*n_i)];
A(i*n_i,:)=1/m*[zeros(1,n_i*(i-2)),k,d,-2*k,-2*d,k,d,zeros(1,n_x-n_i*(i+1))];
end
%i=M, measured output
A(N*n_i-1,:)=[zeros(1,N*n_i-1),1];
A(N*n_i,:)=1/m*[zeros(1,n_x-2*n_i),k,d,-k,-d]; 
%acutation on last mass
n_u=1;
B=1/m*[zeros(n_x-n_i,1);0;1;];
%output: position last mass
C=[1,0,zeros(1,n_i*(N-1))]; 
n_y=size(C,1);
%disturbance: every velocity
E=[];%kron(eye(N),[0;1]);
n_w=size(E,2);
sys = ss(A,[B,E],C,0);