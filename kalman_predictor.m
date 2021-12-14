% MATLAB program for implementing kalman predictor
function main
clear all;
close all;
% initialize system and simulation parameters
global A B   
N=25;
Qf=[1 0;0 1];
Q=[1 0;0 1]; J0=0;
R=1;s=1; n=2;I=eye(2);
x0 = [10+2.5*randn(1);5+2.5*randn(1)]; xk0=[10;5];
A=[0.5 0;-1 1.5];
B=[0.5;0.1];
C=[1 0.5];

% Initializing the vectors and matrices
xk = zeros(2,N+1);
xk(:,1)=xk0;
P=[1 0;0 1];
Pk=zeros(2,2*(N+1));
P0=zeros(2,N+1);
P0(1,1)=P(1,1);
P0(2,1)=P(2,2);
Pk(:,1:2)=P;
Lk=zeros(2,N);
x = zeros(2,N+1);
x(:,1)=x0;
u = zeros(1,N);
y = zeros(1,N);
    

% kalman predictor algorithm      
for j=1:N   
    d(:,j)=0.25*randn(2,1);
    v(j)=0.25*randn(1);
    %K=K0(j,:)
    K=[2.735 -2.747];
    u(j)=-K*x(:,j);
    x(:,j+1)=A*x(:,j)+B*u(1,j)+d(:,j);
    y(j)=C*x(:,j)+v(j);
    Lk(:,j)=A*Pk(:,2*(j)-1:2*(j))*C'*(C*Pk(:,2*(j)-1:2*(j))*C'+R)^-1;  
    xk(:,j+1)= A*xk(:,j)+B*u(j)+Lk(:,j)*(y(j)-C*xk(:,j));
    Pk(:,2*(j+1)-1:2*(j+1))=(A-Lk(:,j)*C)*Pk(:,2*j-1:2*j)*(A-Lk(:,j)*C)'+Lk(:,j)*R*Lk(:,j)'+Q;
    P0(1,j+1)=Pk(1,2*(j+1)-1);
    P0(2,j+1)=Pk(2,2*(j+1));
end      
    
% Plotting the responces  
time = (0:N);
subplot(3,1,1)
plot(time,x(1,:),'k.-',time,x(2,:),'r.-','LineWidth',1) 
hold on
plot(time,xk(1,:),'k.-',time,xk(2,:),'r.-','LineWidth',1) 
legend('$x_{1},\hat{x}_{1}$','$x_{2},\hat{x}_{2}$','Interpreter','latex');
%axis([0 50 -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{x}_{k},\hat{\textbf{x}}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(3,1,2)
plot(time(1:end-1),Lk(1,:),'k.-','LineWidth',1) 
hold on
plot(time(1:end-1),Lk(2,:),'r.-','LineWidth',1) 
legend('$L_{1}$','$L_{2}$','Interpreter','latex');
%axis([0 50 -3 3])
xlabel('k','Interpreter','latex');ylabel('$\textbf{L}_{k}$','Interpreter','latex');
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
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg fig1
end

