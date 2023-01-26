% -- Kuramoto network
% Each node has a phase, theta, following
% thetai = omegai + sum_j=1^n Aij sin(thetai-thetaj)

% supporting files:
% - sparsify dynamics.m
% - pool data.m
% - sparseGalerkin.m
% - kura_rhs.m
clear all; close all; clc
%% Plot behavior of one oscillator
t = linspace(0,200,10000); w = 0.5;
theta = w*t;

plot(t, cos(theta),'k','LineWidth',2)
axis off

%% Construct Kuramoto System
t=linspace(0,200,10000);
dt=t(2)-t(1);
rng(447);
K=10;%10;  % coupling strength
n=100; % number of oscillators
rad=ones(n,1);
thetai=3.14*randn(n,1)%0.4*(2*randn(n,1)); % initial condition
omega=0.4*(rand(n,1)+0.5); % random natural frequency for each oscillator

A=rand(n,n);  
A=(A>0.8).*A; % connectivity matrix
%% Simulate
opts = odeset('RelTol',1e-8,'AbsTol',1e-5);

[tsys,y]=ode45('kura_rhs',t,thetai,opts,omega,n,K,A);

%-- Visualize everything at once
figure(1)
h = pcolor(cos(y).');
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
xticklabels('')
set(gca, 'FontSize', 18,'FontName', 'Cambria')
%% SVD of trajectory
X = cos(y).';
start = 1000;
[U, S, V] = svd(X(:,start:end), 'econ');
% -- Low dimensional dynamics
figure
hold on
plot3(U(:,1).'*X,U(:,2).'*X,U(:,3).'*X,'k','LineWidth',1.5);
plot3(U(:,1).'*X(:,1),U(:,2).'*X(:,1),U(:,3).'*X(:,1),'k.','MarkerSize',10);
set(gca, 'FontSize', 18,'FontName', 'Cambria')
grid on

% -- Energy captured in first 10 modes
figure
plot(diag(S(1:10,1:10))./sum(diag(S(1:10,1:10))),'ro','LineWidth',1.5, 'MarkerFaceColor','w');
set(gca, 'FontSize', 18,'FontName', 'Cambria')

%% SINDy
x1 = U(:,1).'*X(:,start:end); 
x2 = U(:,2).'*X(:,start:end);
x3 = U(:,3).'*X(:,start:end);

% -- calculate derivatives
x1dot = (x1(3:end)-x1(1:end-2))/(2*dt);
x2dot = (x2(3:end)-x2(1:end-2))/(2*dt);

% -- build library from center points of X
A = poolData([x1(2:end-1).',x2(2:end-1).'],2,5,0);

% -- STLSQ to get coefficients
coeffs_x1 = sparsifyDynamics(A,x1dot.',10^(-3),1); 
coeffs_x2 = sparsifyDynamics(A,x2dot.',10^(-3),1);

figure
subplot(1,2,1)
bar(coeffs_x1)
subplot(1,2,2)
bar(coeffs_x2)

Xi = [coeffs_x1,coeffs_x2];

%% Simulate Dynamics
init = [U(:,1).'*X(:,start); U(:,2).'*X(:,start)];

tspan = [0 200]; 
[t,y] = ode45('sparseGalerkin',tspan,init,[],Xi,5,0);

colors = lines(6);
figure
subplot(4,1,1)
hold on
plot(tsys(start:end)-tsys(start),U(:,1).'*X(:,start:end),'color', [colors(1,:), 0.4],'LineWidth',1.8)
plot(t,y(:,1),':', 'color', colors(1,:),'LineWidth',2)
subplot(4,1,2)
hold on
plot(tsys(start:end)-tsys(start),U(:,2).'*X(:,start:end),'color', [colors(2,:), 0.4],'LineWidth',1.8)
plot(t,y(:,2),':', 'color', colors(2,:),'LineWidth',2)
subplot(4,1,[3 4])
hold on
plot(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end),'color', [0 0 0 0.4],'LineWidth',1.8)
plot(y(:,1), y(:,2),':','color',[0 0 0], 'LineWidth',2)
