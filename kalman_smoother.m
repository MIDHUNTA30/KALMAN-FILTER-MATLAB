% MATLAB program for implementing kalman smoother
function main
clear all;
close all;
global A B
% initialize system and simulation parameters
N=25;
Qf=[1 0;0 1];
Q=[1 0;0 1]; J0=0;
R=1;s=1; n=2;I=eye(2);
x0 = [10+2.5*randn(1);5+2.5*randn(1)]; xk0=[10;5];
A=[0.5 0;-1 1.5];
B=[0.5;0.1];
C=[1 0.5];

% Initializing the vectors and matrices
xp = zeros(2,N);
xk = zeros(2,N+1);
xN = zeros(2,N+1);
xk(:,1)=xk0;
P=[1 0;0 1];
Pp=zeros(2,2*(N)); 
Pk=zeros(2,2*(N+1));
PN=zeros(2,2*(N+1));
Pk(:,1:2)=P;
P0=zeros(2,N+1);
L0=zeros(2,N);
Lk=zeros(2,N);
Lsk=zeros(2,N);  

x = zeros(2,N+1);
x(:,1)=x0;
u = zeros(1,N);
y = zeros(1,N);
    
% Kalman smoother algorithm stage 1: Filtering
for j=2:N+1   
    d(:,j-1)=0.25*randn(2,1);
    v(j-1)=0.25*randn(1);
    
    %K=K0(j,:)
    K=[2.735 -2.747];
    u(j-1)=-K*x(:,j-1);
    x(:,j)=A*x(:,j-1)+B*u(1,j-1)+d(:,j-1);
    y(j-1)=C*x(:,j-1)+v(j-1);
      
    xp(:,j-1)=A*xk(:,j-1)+B*u(1,j-1);
    Pp(:,2*(j-1)-1:2*(j-1))=A*Pk(:,2*(j-1)-1:2*(j-1))*A'+Q;

    Lk(:,j-1)=Pp(:,2*(j-1)-1:2*(j-1))*C'*(C*Pp(:,2*(j-1)-1:2*(j-1))*C'+R)^-1;   
    xk(:,j)= xp(:,j-1)+Lk(:,j-1)*(y(j-1)-C*xp(:,j-1));
    Pk(:,2*j-1:2*j)=(I-Lk(:,j-1)*C)*Pp(:,2*(j-1)-1:2*(j-1))*(I-Lk(:,j-1)*C)'+Lk(:,j-1)*R*Lk(:,j-1)';
end      
    
   xN(:,N+1)=xk(:,N+1);
   PN(:,2*(N+1)-1:2*(N+1))=Pk(:,2*(N+1)-1:2*(N+1));
   P0(1,N+1)=PN(1,2*(N+1)-1);
   P0(2,N+1)=PN(2,2*(N+1)); 

% Kalman smoother algorithm stage 2: Smoothing   
for k=N:-1:1
    LN(:,2*k-1:2*k)= Pk(:,2*k-1:2*k)*A'*(A*Pk(:,2*k-1:2*k)*A'+Q)^-1;   
    xN(:,k)=xk(:,k)+LN(:,2*k-1:2*k)*(xN(:,k+1)-xp(:,k));    
    PN(:,2*k-1:2*k)=Pk(:,2*k-1:2*k)+LN(:,2*k-1:2*k)*(PN(:,2*(k+1)-1:2*(k+1))-Pp(:,2*k-1:2*k))*LN(:,2*k-1:2*k)';
    P0(1,k)=PN(1,2*k-1);
    P0(2,k)=PN(2,2*k); 
    L0(1,k)=LN(1,2*k-1);
    L0(2,k)=LN(2,2*k); 
end    

% Plotting the responces    
figure(1)
time = (0:N);
subplot(3,1,1)
plot(time,x(1,:),'k.-',time,x(2,:),'r.-','LineWidth',1) 
hold on
plot(time,xN(1,:),'k.-',time,xN(2,:),'r.-','LineWidth',1) 
legend('$x_{1},\hat{x}_{1}$','$x_{2},\hat{x}_{2}$','Interpreter','latex');
%axis([0 50 -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{x}_{k},\hat{\textbf{x}}_{k|N}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(3,1,2)
plot(time(1:end-1),L0(1,:),'k.-','LineWidth',1) 
hold on
plot(time(1:end-1),L0(2,:),'r.-','LineWidth',1) 
legend('$L_{s_{11}}$','$L_{s_{22}}$','Interpreter','latex');
%axis([0 50 -3 3])
xlabel('k','Interpreter','latex');ylabel('$\textbf{L}_{s_{k}}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(3,1,3)
plot(time,P0(1,:),'k.-','LineWidth',1) 
hold on
plot(time,P0(2,:),'r.-','LineWidth',1)
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
%axis([0 50 -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k|N}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg fig3
end

