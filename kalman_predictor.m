% MATLAB program for implementing kalman predictor
clear all;
close all;
% initialize system and simulation parameters
A=[0.5 0;-1 1.5];
B=[0.5;0.1];
C=[1 0.5];
N=25;n=2;m=1;p=1;I=eye(n); 
Q=[1 0;0 1];R=1;
K=[2.735 -2.747]; 


% Initializing the vectors and matrices
x0 = [10+2.5*randn(1);5+2.5*randn(1)]; x0a=[10;5];
xp = zeros(n,N+1); % xp stores the sequence of predicted states 
xp(:,1)=x0a;
x = zeros(n,N+1);  % x stores the sequence of actual states xk
x(:,1)=x0;
u = zeros(m,N);
y = zeros(p,N);
Pp=zeros(n,n*(N+1));  % Pp stores the sequence of matrices Pk
Pp(:,1:2)=Q; 
Ppd=zeros(n,N+1);  % Ppd stores the diagonal elements of Pk
Ppd(:,1)=diag(Q);
L=zeros(n,N);  % L stores Lk sequence
  

% kalman predictor algorithm      
for j=1:N   
    d(:,j)=0.25*randn(n,1);
    v(:,j)=0.25*randn(p);
    u(:,j)=-K*x(:,j);  % Here the actual state is used in the control law since the system is open-loop unstable. In other examples the estimated state can be used.
    x(:,j+1)=A*x(:,j)+B*u(:,j)+d(:,j);
    y(:,j)=C*x(:,j)+v(:,j);
    L(:,j)=A*Pp(:,n*(j-1)+1:n*j)*C'*(C*Pp(:,n*(j-1)+1:n*j)*C'+R)^-1;  
    xp(:,j+1)= A*xp(:,j)+B*u(j)+L(:,j)*(y(j)-C*xp(:,j));
    Pp(:,n*j+1:n*(j+1))=(A-L(:,j)*C)*Pp(:,n*(j-1)+1:n*j)*(A-L(:,j)*C)'+L(:,j)*R*L(:,j)'+Q;
    Ppd(:,j+1)=diag(Pp(:,n*j+1:n*(j+1)));
end      
    
% Plotting the responces  
time = (0:N);
subplot(3,1,1)
plot(time,x(1,:),'k.-',time,x(2,:),'r.-','LineWidth',1) 
hold on
plot(time,xp(1,:),'k.-',time,xp(2,:),'r.-','LineWidth',1) 
legend('$x_{1},\hat{x}_{1}$','$x_{2},\hat{x}_{2}$','Interpreter','latex');
%axis([0 50 -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{x}_{k},\hat{\textbf{x}}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(3,1,2)
plot(time(1:end-1),L(1,:),'k.-',time(1:end-1),L(2,:),'r.-','LineWidth',1) 
legend('$L_{1}$','$L_{2}$','Interpreter','latex');
%axis([0 50 -3 3])
xlabel('k','Interpreter','latex');ylabel('$\textbf{L}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(3,1,3)
plot(time,Ppd(1,:),'k.-',time,Ppd(2,:),'r.-','LineWidth',1) 
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
%axis([0 50 -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg fig1

