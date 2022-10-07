% MATLAB program for implementing kalman smoother
clear all;
close all;
% initialize system and simulation parameters
A=[0.5 0;-1 1.5];
B=[0.5;0.1];
C=[1 0.5];
N=20; n=2;m=1;p=1;I=eye(n); 
Q=[1 0;0 1];R=1;
K=[2.735 -2.747];


% Initializing the vectors and matrices
x0 = [10+2.5*randn(1);5+2.5*randn(1)]; x0a=[10;5];
x = zeros(n,N+1);  % x stores the sequence of actual states xk
x(:,1)=x0;
u = zeros(m,N);
y = zeros(p,N);
xp = zeros(n,N);   % xp stores the sequence of predicted states xk|k-1
xf = zeros(n,N+1);   % xf stores the sequence of filtered states xk|k
xf(:,1)=x0a;
xs = zeros(n,N+1);   % xs stores the sequence of smoothed states xk|N
Pp=zeros(n,n*(N));  % Pp stores the sequence of matrices Pk|k-1
Pf=zeros(n,n*(N+1));  % Pf stores the sequence of matrices Pk|k
Pf(:,1:n)=Q;
Ps=zeros(n,n*(N+1));  % Ps stores the sequence of matrices Pk|N
Psd=zeros(n,N+1);  % Psd stores the diagonal elements of Pk|N
L=zeros(n,N);  % L stores the sequence of gains Lk
Ls=zeros(n,n*N);  % Ls stores the sequence of gains Lsk
Lsd=zeros(n,N);  % Lsd stores the diagonal elements of Lsk


    
% Kalman smoother algorithm stage 1: Filtering
for j=2:N+1   
    d(:,j-1)=0.25*randn(n,1);
    v(j-1)=0.25*randn(p);
    u(j-1)=-K*x(:,j-1);
    Ak=A+((-1)^j)*0.5*eye(n);Bk=B+((-1)^j)*0.1*B;
    x(:,j)=Ak*x(:,j-1)+Bk*u(1,j-1)+d(:,j-1);
    y(j-1)=C*x(:,j-1)+v(j-1);
    xp(:,j-1)=Ak*xf(:,j-1)+Bk*u(1,j-1);
    Pp(:,n*(j-2)+1:n*(j-1))=Ak*Pf(:,n*(j-2)+1:n*(j-1))*Ak'+Q;
    L(:,j-1)=Pp(:,n*(j-2)+1:n*(j-1))*C'*(C*Pp(:,n*(j-2)+1:n*(j-1))*C'+R)^-1;   
    xf(:,j)= xp(:,j-1)+L(:,j-1)*(y(j-1)-C*xp(:,j-1));
    Pf(:,n*(j-1)+1:n*j)=(I-L(:,j-1)*C)*Pp(:,n*(j-2)+1:n*(j-1))*(I-L(:,j-1)*C)'+L(:,j-1)*R*L(:,j-1)';
end      
    
   xs(:,N+1)=xf(:,N+1);
   Ps(:,n*N+1:n*(N+1))=Pf(:,n*N+1:n*(N+1));
   Psd(:,N+1)=diag(Ps(:,n*N+1:n*(N+1)));

% Kalman smoother algorithm stage 2: Smoothing   
for k=N:-1:1
    Ls(:,n*(k-1)+1:n*k)= Pf(:,n*(k-1)+1:n*k)*A'*(A*Pf(:,n*(k-1)+1:n*k)*A'+Q)^-1;   
    xs(:,k)=xf(:,k)+Ls(:,n*(k-1)+1:n*k)*(xs(:,k+1)-xp(:,k));    
    Ps(:,n*(k-1)+1:n*k)=Pf(:,n*(k-1)+1:n*k)+Ls(:,n*(k-1)+1:n*k)*(Ps(:,n*k+1:n*(k+1))-Pp(:,n*(k-1)+1:n*k))*Ls(:,n*(k-1)+1:n*k)';
    Psd(:,k)=diag(Ps(:,n*(k-1)+1:n*k));
    Lsd(:,k)=diag(Ls(:,n*(k-1)+1:n*k)); 
end    

% Plotting the responces  
time = (0:N);
subplot(2,2,1)
plot(time,x(1,:),'k.-',time,xs(1,:),'r.-','LineWidth',1) 
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
plot(time,x(2,:),'k.-',time,xs(2,:),'r.-','LineWidth',1) 
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
plot(time(1:end-1),Lsd(1,:),'k.-',time(1:end-1),Lsd(2,:),'r.-','LineWidth',1) 
legend('$L_{1}$','$L_{2}$','Interpreter','latex');
axis([0 20 0 1])
xlabel('k','Interpreter','latex');ylabel('$\textbf{L}_{s_k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
set(gca,'xtick',[0:5:20])
set(gca,'ytick',[0:.25:1])
subplot(2,2,4)
plot(time,Psd(1,:),'k.-',time,Psd(2,:),'r.-','LineWidth',1) 
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
axis([0 20 0 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k|N}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
set(gca,'xtick',[0:5:20])
set(gca,'ytick',[0:2.5:10])
print -dsvg fig2