% -- Figure 7: 4 Examples
% Heterogeneous Oscillators
clear all; close all; clc
%% Set up Parameters

% -- some examples with cool results
%rng(72) %with nK 30, (SVD) start = 10000, (fig) start = 13000
%rng(1024) %with nK 53, (SVD) start = 10000, (fig) start = 1
rng(30) %with nK 65, (SVD) start = 10000, (fig) start = 10000
%rng(200) %with nK 80, (SVD) start = 10000, (fig) start = 5000

n = 100; % number of total oscillators
A=rand(n,n);  
A=(A>0.8).*A; % connectivity

% -- Kuramoto parameters
a =.6;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 10;
nK=65; % number of Kuramoto oscillators
omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each oscillator
% -- FHN parameters 
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;

%% Simulate
 %-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];
% -- time domain
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);
% -- simulate
opts = odeset('RelTol',1e-8,'AbsTol',1e-5);
[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);

% -- raster plot full time series
figure
h = pcolor([cos(y(end-5000:end,1:nK)) y(end-5000:end,nK+1:n)].');
xticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')

% -- raster plot final quarter
figure
h = pcolor([cos(y(:,1:nK)) y(:,nK+1:n)].');
xticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')

%% SVD
start = 10000;
X = [cos(y(:,1:nK)) y(:,nK+1:n) y(:,n+1:end)].';
[U1,S1,V1] = svd(X(:,start:end),'econ');
temp = find(sqrt(cumsum(diag(S1).^2))./sqrt(sum(diag(S1).^2))>.95);
dim = temp(1)
  
% Stagger start further if want to avoid clutter of transient initial dynamics
start = 10000; 

figure
hold on
plot3(U1(:,1).'*X(:,start:end),U1(:,2).'*X(:,start:end), U1(:,3).'*X(:,start:end),'k', 'LineWidth',1.5)
%plot3(U1(:,1).'*X(:,start),U1(:,2).'*X(:,start), U1(:,3).'*X(:,start),'g.', 'MarkerSize',10)
%plot3(U1(:,1).'*X(:,end),U1(:,2).'*X(:,end), U1(:,3).'*X(:,end),'r.', 'MarkerSize',10)
set(gca, 'FontSize',20, 'FontName', 'Cambria')