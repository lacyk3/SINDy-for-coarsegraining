% Figure 8: SINDy for coarsegraining heterogeneous networks
% Apply SINDy to two examples of mixed Kuramoto-FHN networks

% -- supporting files
% - poolData.m
% - SR3.m
% - sparseGalerkin.m


clear all; close all; clc
%% Example (a)
rng(200) 
n = 100; % number of total oscillators
A=rand(n,n);  
A=(A>0.8).*A; % connectivity

% -- Kuramoto parameters
a =.6;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 10;
nK=80; % number of Kuramoto oscillators
omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each oscillator
% -- FHN parameters 
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;

%%-- Simulate
 %-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);
% -- Simulate until system is on limit cycle
opts = odeset('RelTol',1e-8,'AbsTol',1e-5);
[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);
% --- Simulate at high resolution from end of previous sim
tspan =linspace(0,140,2^17);
dt = tspan(2) -tspan(1);
x = y(end,:).';
opts = odeset('RelTol',1e-8,'AbsTol',1e-5);
[t1,y1] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);

%%-- SVD + differentiate
X1 = [cos(y1(:,1:nK)) y1(:,nK+1:end)].'; 
[U1,S1,V1] = svd(X1,'econ');
X = [(U1(:,1).'*X1).', (U1(:,2).'*X1).',(U1(:,3).'*X1).'];
Xs = X(2:end-1,:);
dxdt= (X(3:end,:)-X(1:end-2,:))/(2*dt);
ts = t1(2:end-1);

%%-- SINDy
LIB = poolData(Xs,3,2,0); 
kappa = 80; 
tau = 10^(-8);
ss= .01; %step size
tf = 0;%trimmed fraction
maxiter = 600;

XiF = []; XiS = [];
for j = 1:3
[Xi_full,Xi_sparse,trimmed_points] = SR3(LIB, dxdt(:,j),tau, kappa, ss, tf, maxiter);
XiF(:,j) = Xi_full;
XiS(:,j) = Xi_sparse;
end

%-- Simulate
init = Xs(1,:).';
tspan = [0 100]; 
[ty,y] = ode45('sparseGalerkin',tspan,init,opts,XiS,2,0);

%--Plot
c = lines
figure
subplot(3,1,1)
hold on
plot(ts,Xs(:,1), 'color',[c(1,:) 0.3], 'LineWidth',1.6)
plot(ty,y(:,1),':', 'color',c(1,:), 'LineWidth',2)
subplot(3,1,2)
hold on
plot(ts,Xs(:,2), 'color',[c(2,:) 0.3], 'LineWidth',1.6)
plot(ty,y(:,2),':', 'color',c(2,:), 'LineWidth',2)
subplot(3,1,3)
hold on
plot(ts,Xs(:,3), 'color',[c(5,:) 0.3], 'LineWidth',1.6)
plot(ty,y(:,3),':', 'color',c(5,:), 'LineWidth',2)
%% Example(b)
rng(2017)
n = 100; % number of total oscillators
A=rand(n,n);  
A=(A>0.8).*A; % connectivity

% -- Kuramoto parameters
a =.6;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 10;
nK=40; % number of Kuramoto oscillators
omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each oscillator
% -- FHN parameters 
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;

%%-- Simulate
 %-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);
% -- Simulate until system is on limit cycle
opts = odeset('RelTol',1e-8,'AbsTol',1e-5);
[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);
% --- Simulate at high resolution from end of previous sim
tspan =linspace(0,350,2^17);
dt = tspan(2) -tspan(1);
x = y(end,:).';
opts = odeset('RelTol',1e-8,'AbsTol',1e-5);
[t1,y1] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);
X1 = [cos(y1(:,1:nK)) y1(:,nK+1:end)].';

%%
[U1,S1,V1] = svd(X1(n:end,:),'econ');
[U2,S2,V2] = svd(X1(nK:n,:),'econ');
[U3,S3,V3] = svd(X1(1:nK,:),'econ');

Xs = [(U1(:,1).'*X1(n:end,:)).', (U1(:,2).'*X1(n:end,:)).',(U1(:,3).'*X1(n:end,:)).',...
    (U2(:,1).'*X1(nK:n,:)).', (U2(:,2).'*X1(nK:n,:)).',(U2(:,3).'*X1(nK:n,:)).',...
    (U3(:,1).'*X1(1:nK,:)).', (U3(:,2).'*X1(1:nK,:)).',(U3(:,3).'*X1(1:nK,:)).'];
dxdt= (Xs(3:end,:)-Xs(1:end-2,:))/(2*dt);

figure
plot3(Xs(:,1),Xs(:,2),Xs(:,3), 'k','LineWidth', 1.8)
grid on
view([-18.1333   16.9333])
set(gca, 'FontSize', 18, 'FontName', 'Cambria')
figure
plot3(Xs(:,4),Xs(:,5),Xs(:,6),'k','LineWidth', 1.8)
grid on
set(gca, 'FontSize', 18, 'FontName', 'Cambria')
view([0.5333   90.0000])
figure
plot3(Xs(:,7),Xs(:,8),Xs(:,9),'k','LineWidth', 1.8)
view([57.8666   19.6000])
grid on
set(gca, 'FontSize', 18, 'FontName', 'Cambria')
%%

temp1 = poolData(Xs(2:end-1,4:6),3,3,0);
temp2 = poolData(Xs(2:end-1,7:9),3,3,0);
LIB = [poolData(Xs(2:end-1,1:3),3,3,0) temp1(:,2:end) temp2(:,2:end)];
kappa = 80; 
tau = 10^(-8);
ss= .01; %step size
tf = 0;%trimmed fraction
maxiter = 600;
% 
XiF = []; XiS = [];
for j = 1:9
[Xi_full,Xi_sparse,trimmed_points] = SR3(LIB, dxdt(:,j),tau, kappa, ss, tf, maxiter);
XiF(:,j) = Xi_full;
XiS(:,j) = Xi_sparse;
end
% 
figure
colors = lines;
for j = 1:9
subplot(3,3,j)
hold on
plot(t1(2:end-1),dxdt(:,j), 'color', [colors(j,:) 0.3], 'LineWidth',1.5)
plot(t1(2:end-1),LIB*XiF(:,j),':', 'color', colors(j,:) , 'LineWidth',2)
end

figure
imagesc([abs(XiS(:,1:3))*10 abs(XiS(:,4:6)) abs(XiS(:,7:9))*10])
caxis([0 100])
c = gray;
c = flipud(c);
colormap(c);
axis off