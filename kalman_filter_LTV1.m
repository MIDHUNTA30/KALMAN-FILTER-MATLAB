% MATLAB program for implementing kalman filter
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
x = zeros(n,N+1);  % x stores the sequence of actual states xk
x(:,1)=x0;
u = zeros(m,N);
y = zeros(p,N);
xp = zeros(n,N);  % xp stores the sequence of predicted states xk|k-1
xf = zeros(n,N+1);  % xf stores the sequence of filtered states xk|k
xf(:,1)=x0a;
Pp=zeros(n,n*(N));  % Pp stores the sequence of matrices Pk|k-1
Pf=zeros(n,n*(N+1));  % Pf stores the sequence of matrices Pk|k
Pf(:,1:2)=Q;
Pfd=zeros(n,N+1); % Pfd stores the diagonal elements of Pk|k
Pfd(:,1)=diag(Q);
L=zeros(n,N);  % L stores the sequence of gains Lk
    

% kalman filter algorithm          
for j=2:N+1   
    d(:,j-1)=0.25*randn(n,1);
    v(:,j-1)=0.25*randn(p);
    u(j-1)=-K*x(:,j-1);
    Ak=A+((-1)^j)*0.5*eye(n);Bk=B+((-1)^j)*0.1*B;
    x(:,j)=Ak*x(:,j-1)+Bk*u(:,j-1)+d(:,j-1);
    y(:,j-1)=C*x(:,j-1)+v(:,j-1);
    xp(:,j-1)=Ak*xf(:,j-1)+Bk*u(:,j-1);
    Pp(:,n*(j-2)+1:n*(j-1))=Ak*Pf(:,n*(j-2)+1:n*(j-1))*Ak'+Q;
    L(:,j-1)=Pp(:,n*(j-2)+1:n*(j-1))*C'*(C*Pp(:,n*(j-2)+1:n*(j-1))*C'+R)^-1;  
    xf(:,j)= xp(:,j-1)+L(:,j-1)*(y(j-1)-C*xp(:,j-1));
    Pf(:,n*(j-1)+1:n*j)=(I-L(:,j-1)*C)*Pp(:,n*(j-2)+1:n*(j-1))*(I-L(:,j-1)*C)'+L(:,j-1)*R*L(:,j-1)';
    Pfd(:,j)=diag(Pf(:,n*(j-1)+1:n*j));
    %u(j)=-K*xf(:,j);
end      
    
% Plotting the responces  
time = (0:N);
subplot(2,2,1)
plot(time,x(1,:),'k.-',time,xf(1,:),'r.-','LineWidth',1) 
legend('$x_{1}$','$\hat{x}_{1}$','Interpreter','latex');
axis([0 20 -10 10])
xlabel('k','Interpreter','latex');ylabel('${x}_{1_k},\hat{{x}}_{1_k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
set(gca,'xtick',[0:5:20])
set(gca,'ytick',[-10:5:10])
subplot(2,2,2)
plot(time,x(2,:),'k.-',time,xf(2,:),'r.-','LineWidth',1) 
legend('$x_{2}$','$\hat{x}_{2}$','Interpreter','latex');
axis([0 20 -10 10])
xlabel('k','Interpreter','latex');ylabel('${x}_{2_k},\hat{{x}}_{2_k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
set(gca,'xtick',[0:5:20])
set(gca,'ytick',[-10:5:10])
subplot(2,2,3)
plot(time(1:end-1),L(1,:),'k.-',time(1:end-1),L(2,:),'r.-','LineWidth',1) 
legend('$L_{1}$','$L_{2}$','Interpreter','latex');
axis([0 20 0 2])
xlabel('k','Interpreter','latex');ylabel('$\textbf{L}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
set(gca,'xtick',[0:5:20])
set(gca,'ytick',[0:.5:2])
subplot(2,2,4)
plot(time,Pfd(1,:),'k.-',time,Pfd(2,:),'r.-','LineWidth',1) 
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
axis([0 20 0 20])
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k|k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
set(gca,'xtick',[0:5:20])
set(gca,'ytick',[0:5:20])
print -dsvg fig2
