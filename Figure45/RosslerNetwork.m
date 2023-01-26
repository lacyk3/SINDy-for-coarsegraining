% -- Rossler network
% The state of each node is defined by three state variables, x, y, and z,
% which follow

% xi = - yi - zi + (K/n)sum_j=1^n Aijsin(xj-xi)
% yi = xi + ayi + (K/n)sum_j=1^n Aijsin(yj-yi)
% zi = b + zi(xi-c) + (K/n)sum_j=1^n Aijsin(zj-zi)

% -- supporting files
% - Rossler_rhs.m
% - poolData.m
% - SR3.m
% - sparseGalerkin.m

clear all; close all; clc
%% Plot behavior of one oscillator
% a = 0.2; b = 0.2; c = 25;
% 
% Rossler = @(t,x)([ -x(2) - x(3); 
%                     x(1) + a*x(2);
%                     b + x(3).*(x(1)-c)]);
%                 
% options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
% [t,y] = ode45(Rossler,[0 100],[0;1;1],options);
% 
% figure
% hold on
% plot(t, y(:,1),'k-.', 'LineWidth',1.5)
% plot(t, y(:,2), 'k:', 'LineWidth', 1.5)
% plot(t, y(:,3), 'k', 'LineWidth', 1.5)
% hold off
% axis([0 100 -75 250])
% set(gca, 'FontSize', 18,'FontName', 'Cambria')
%% Simulate network
rng(1)
n=50;
alpha = [0.1, 0.1, 14];
%A = ones(n,n);
A=rand(n,n);  
A=(A>0.8).*A; 
K = 90;
on = 1;

x0 = 2*randn(3*n,1);                
options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
tspan = linspace(0,500,2^13);
dt = tspan(2) - tspan(1);
[t1,x1] = ode45('Rossler_rhs',tspan,x0,options,alpha,n,K,A,on);
tspan = linspace(0,100,2^14);
dt = tspan(2) - tspan(1);
[t,x] = ode45('Rossler_rhs',tspan,x1(end,:).',options,alpha,n,K,A,on);

figure
h = pcolor(x.'/10);
set(h,'EdgeColor','none');
colormap('jet')
colorbar()
xticklabels('')
set(gca, 'FontSize',18,'FontName','Cambria')
%% SVD
% The dynamics live on a tilted version of the attractor because everything
% syncs up

start = 1;
[U,S,V] = svd(x(start:end,:).', 'econ');

figure
plot(diag(S(1:10,1:10))./cumsum(diag(S(1:10,1:10))),'ro','LineWidth',1.5)
set(gca, 'FontSize',18,'FontName','Cambria')
figure
plot3(U(:,1).'*(x.'),U(:,2).'*(x.'),U(:,3).'*(x.'),'k','LineWidth',2)
axis([-200 150 -150 150 -50 250])
grid on
set(gca, 'FontSize',18,'FontName','Cambria')
%% SINDy to discover dynamics
start = 1;
X = [(U(:,1).'*(x(start:end,:).')).',(U(:,2).'*(x(start:end,:).')).', (U(:,3).'*(x(start:end,:).')).'];
Xdot = (X(3:end,:) - X(1:end-2,:))/(2*dt);

x1 = X(:,1); x2 = X(:,2); x3= X(:,3);
x1dot = Xdot(:,1); x2dot = Xdot(:,2); x3dot= Xdot(:,3);

A = poolData(X(2:end-1,:),3,2,0);

kappa = 80; 
tau = 10^(-1);%threshold %10^(-4)
ss= .01; %step size
tf = 0.2;%trimmed fraction
maxiter = 600;


[coeffs_x1, Xi1, trimmed_points1] = SR3(A,x1dot,tau,kappa, ss, tf, maxiter); 
[coeffs_x2, Xi2, trimmed_points2]  = SR3(A,x2dot,tau,kappa, ss, tf, maxiter);
[coeffs_x3, Xi3, trimmed_points3]  = SR3(A,x3dot,tau,kappa, ss, tf, maxiter); 

figure
hold on
plot3(x1(2:end-1),x2(2:end-1),x3(2:end-1),'k', 'LineWidth', 2)
plot3(x1(trimmed_points3<0.8015),x2(trimmed_points3<0.8015),x3(trimmed_points3<0.8015),'r.','MarkerSize',15)
axis([-200 150 -150 150 -50 250])
grid on
set(gca, 'FontSize',18,'FontName','Cambria')
view([143.4667   28.1333])

figure
subplot(1,3,1)
bar(coeffs_x1)
subplot(1,3,2)
bar(coeffs_x2)
subplot(1,3,3)
bar(coeffs_x3)

figure
subplot(3,1,1)
plot(Xdot(:,1))
hold on
plot(A*coeffs_x1)
subplot(3,1,2)
plot(Xdot(:,2))
hold on
plot(A*coeffs_x2)
subplot(3,1,3)
plot(Xdot(:,3))
hold on
plot(A*coeffs_x3)

Xi = zeros(numel(coeffs_x3),3);
Xi(:,1) = coeffs_x1; 
Xi(:,2) = coeffs_x2; 
Xi(:,3) = coeffs_x3;


figure
imagesc(abs([[Xi(1:4,1),Xi(1:4,2), Xi(1:4,3)];60*[Xi(5:end,1),Xi(5:end,2), Xi(5:end,3)]] ))
c = gray;
c = flipud(c);
colormap(c);
axis off
%% Simulate discovered dynamics

init = X(1,:).';
tspan = [0 100]; 
[ty,y] = ode45('sparseGalerkin',tspan,init,options,Xi,2,0);

colors = lines(5);

figure
subplot(3,1,1)
hold on
plot(t(1:10:end),X(1:10:end,1), 'color', [colors(1,:), 0.4],'LineWidth',1.8)
plot(ty(1:10:end),y(1:10:end,1),':','color', colors(1,:),'LineWidth',2.2)
axis([0 100 -200 200])
xticklabels('')
yticklabels('')
set(gca, 'FontSize',18,'FontName','Cambria')
subplot(3,1,2)
hold on
plot(t(1:10:end),X(1:10:end,2), 'color', [colors(2,:), 0.4],'LineWidth',1.8)
plot(ty(1:10:end), y(1:10:end,2),':','color', colors(2,:),'LineWidth',2.2)
axis([0 100 -200 200])
xticklabels('')
yticklabels('')
set(gca, 'FontSize',18,'FontName','Cambria')
subplot(3,1,3)
hold on
plot(t(1:10:end),X(1:10:end,3),'color', [colors(5,:), 0.4],'LineWidth',1.8)
plot(ty(1:10:end),y(1:10:end,3),':','color', colors(5,:),'LineWidth',2.2)
axis([0 100 -50 250])
yticklabels('')
set(gca, 'FontSize',18,'FontName','Cambria')
