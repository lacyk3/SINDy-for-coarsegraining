% Figure 1

% Let's see if trimming can pick out short periods of impulsive dynmamics
% on a few example systems!

% Necessary supporting files that must be in path: 
% - poolData.m (for library construction)
% - SR3.m (for running SINDy with trimming on one dimensional system)
% - FHN_rhs.m (function for simulating network of FHN oscillators)

clear all; close all; clc
%% Rayleigh Oscillator Example

% The Rayleigh equation can be formulated: eps*y'' - (1-y'^2)y' + y = 0.     
% Making the substitution z = y' we get a system of two coupled ODEs:
% y' = z
% eps*z' = -y + (1-z^2)z

% -- generate data
eps = .001;
init = [0;-1.72];
Rayleigh = @(t,x) ([x(2);
                    1/eps*(x(2) - x(2)^3/3 - x(1))]);
tspan = linspace(0,2.5,10000);
% use tspan with denser cluster of points in boundary layer so that it
% shows up in plot --> use trimmed fraction = 0.2
% tspan =  [linspace(0,.46,800),linspace(.460001,.48,1000),...
%          linspace(.48001,1.31,1000),linspace(1.31001,1.33,1000),...
%          linspace(1.330001,2.14,1000), linspace(2.14001,2.16,1000),...
%          linspace(2.16001,2.5,1000)];     
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45(Rayleigh,tspan,init,options);

% -- sample derivative
dx_centerdiff = (y(3:end,:) - y(1:end-2,:))./...
    [(t(3:end) - t(1:end-2)) (t(3:end) - t(1:end-2))];

% -- SINDy via SR3 with trimming
% build library 
nvars = 1; polyorder = 3; usesine = 0;
A = poolData(y(2:end-1,2),nvars,polyorder,usesine);
% hyperparameters
trimmed_fraction = .011; %0.2;
kappa = 80; 
stepsize = .01; 
maxiter = 2000;
lambda = .9;
% apply SR3 algorithm
[Xi_full,Xi_sparse,trimmed_points] = ...
    SR3(A, dx_centerdiff(:,2),lambda, kappa, stepsize, trimmed_fraction, maxiter);

threshold = 1 - trimmed_fraction; 
figure 
hold on
% plot axes
plot([-1 1],[0 0],'k'); 
plot([0 0],[-2.5 2.5],'k');
% plot outer solution from approximation theory
z_0 = linspace(-3,3,500);
y_out = z_0.*(1-z_0.^2./3 );
plot(y_out,z_0, 'k:','LineWidth',2); hold on
% plot trajectory
plot(y(:,1),y(:,2),'k.', 'LineWidth',2);
% plot trimmed points in red
plot(y(trimmed_points<threshold,1),y(trimmed_points<threshold,2),'r.', 'MarkerSize',10);
hold off

set(gca, 'FontSize', 18,'FontName', 'Cambria')
axis([-1 1 -2.5 2.5])
set(gcf,'position',[100,10,300,400])

%% Rossler Example
% This well-studied system makes infrequent excursions in the z direction, 
% which perhaps we can trim off.

a = 0.2; b = 0.2; c = 25;

Rossler = @(t,x)([ -x(2) - x(3); 
                    x(1) + a*x(2);
                    b + x(3).*(x(1)-c)]);
                
 
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45(Rossler,[0 100],[0;1;1],options);
dx_centerdiff = (y(3:end,:) - y(1:end-2,:))./...
    [(t(3:end) - t(1:end-2)) (t(3:end) - t(1:end-2)) (t(3:end) - t(1:end-2))];

figure
hold on
plot3(y(:,1),y(:,2),y(:,3),'k','LineWidth',1.8);
set(gca, 'FontSize', 18,'FontName', 'Cambria')
grid on

nvars = 1; polyorder = 2; usesine = 0;
A = poolData(y(2:end-1,3),nvars,polyorder,usesine);

kappa = 80; 
stepsize = .01; 
trimmed_fraction = 0.2;
lambda = 0.8;
maxiter = 2000;
[Xi_full,Xi_sparse,trimmed_points] = ...
    SR3(A, dx_centerdiff(:,3),lambda, kappa, stepsize, trimmed_fraction, maxiter);
threshold =1-trimmed_fraction;
% -- Visualize Results
plot3(y(trimmed_points<threshold,1),...
    y(trimmed_points<threshold,2),...
    y(trimmed_points<threshold,3),'r.', 'MarkerSize',10);
hold off

%% FHN Network Example
% In this networked dynamical system, each node is a simplified model of a 
% neuron firing. For a single oscillator, the neuron will not spike
% continuously without an external stimulus, but coupling them together can
% provide the necessary input such that sustained spiking patterns persist

% -- Construct FHN System
K=.2;  % coupling strength
n=50; % number of oscillators
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
A=ones(n,n);  % all to all connections
Vi= [-.95 + 2*rand(n,1); -.05 + .1*rand(n,1)]; % initial state

% -- simulate "data"
tspan = linspace(0,1000,20000);
opts = odeset('RelTol',10^(-14),'AbsTol',10^(-14));
[t,y]=ode45(@(t,y) FHN_rhs(t,y,[],alpha,n,K,A,1),tspan,Vi);

%-- Visualize voltages 
figure
start = 10000;
h = pcolor(y(start:end,1:n).');
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
xticklabels('')
set(gca, 'FontSize', 18,'FontName', 'Cambria')

% -- SVD
[Uv, Sv, Vv] = svd(y(start:end,1:n).', 'econ');
x1 = Sv(1,1)*Vv(:,1); 
x2 = Sv(2,2)*Vv(:,2); 
x3 = Sv(3,3)*Vv(:,3); 

figure
hold on
plot(x1,x2, 'k', 'LineWidth', 1.5)

% -- derivatives
dt = t1(2) - t1(1);
X = x1;
dXs = (X(3:end,:)-X(1:end-2,:))./(2*dt);
Xs = X(2:end-1,:);

% -- SINDy with trimming
A = poolData(Xs,1,3,0);
kappa = 80; 
stepsize = .01; 
trimmed_fraction = 0.2;
lambda = 0.8;
maxiter = 200;
[Xi_full,Xi_sparse,trimmed_points] = ...
    SR3(A, dXs,lambda, kappa, stepsize, trimmed_fraction, maxiter);
threshold =1-trimmed_fraction;
plot(x1(trimmed_points<threshold),x2(trimmed_points<threshold), 'r.')