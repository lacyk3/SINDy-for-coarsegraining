% -- Fitzhugh Nagumo Network
% Each node has a voltage, v, and a recovery variable, w, following
% vi = a3vi^3 + a2vi^2 a1vi - wi + sum_j=1^j=nAij(vi-vj)
% wi = cvi - bwi

% -- supporting files
% - FHN_rhs.m
% - poolData.m
% - DifferentLibraries_rhs
clear all; close all; clc
%% Plot behavior of one oscillator
I = @(t)(1); % -- include external forcing, otherwise only one spike
a = 0.7; b = 0.8; tau=100;
FHN = @(t,x)([(x(1)-x(1).^3)./3 - x(2) + I(t);
              (x(1) + a - b*x(2))/tau]);

options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45(FHN,[0 400],[-.5;-.5],options);

figure
hold on
plot(t,y(:,1),'k', 'LineWidth',2);
plot(t,y(:,2),'k:', 'LineWidth',2);
set(gca, 'FontSize', 18,'FontName', 'Cambria')
axis off

%% Simulate network
rng(195)
K=.3;  % coupling strength
n=50; % number of oscillators
Vi= [-.95 + 2*rand(n,1); -.05 + .1*rand(n,1)]; % initiate state

% A=rand(n,n);  
% A=(A>0.1).*A;
A=ones(n,n);  % all to all connections
alpha =[-.1,1.1,-1,.01,.01]; 

t0=linspace(0,100,10000);
t1 = linspace(0,2000,40000);

opts = odeset('RelTol',10^(-9),'AbsTol',10^(-9));
[t0,y0]=ode45(@(t,y) FHN_rhs(t,y,opts,alpha,n,K,A,1),t0,Vi);
[t1,y1]=ode45(@(t,y) FHN_rhs(t,y,opts,alpha,n,K,A,1),t1,y0(end,:).');
times= t1;

start = 15000;
%-- Visualize everything at once
figure
h = pcolor(y1(start:end,1:n).');
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
xticklabels('')
set(gca, 'FontSize', 18,'FontName', 'Cambria')
%% SVD
[Uv, Sv, Vv] = svd(y1(start:end,1:n).', 'econ');
[Ur, Sr, Vr] = svd(y1(start:end,n+1:end).', 'econ');

x1 = Sr(1,1)*Vr(:,1); x4 = Sv(1,1)*Vv(:,1);
x2 = Sr(2,2)*Vr(:,2); x5 = Sv(2,2)*Vv(:,2);
x3 = Sr(3,3)*Vr(:,3); x6 = Sv(3,3)*Vv(:,3);

figure
subplot(1,2,1)
plot3(x1,x2,x3,'k','LineWidth',2)
grid on
subplot(1,2,2)
plot3(x4,x5,x6,'k','LineWidth',2)
grid on

%% SINDy to discover dynamics
dt = t1(2) - t1(1);
X = [x1,x2,x3,x4,x5,x6];
dXs = (X(3:end,:)-X(1:end-2,:))./(2*dt);
Xs = X(2:end-1,:);

usesine =0;
polyorder = 1;
nvars = 6;
A = poolData(Xs,nvars,polyorder,usesine);
[m,n] = size(A);
% -- Set regularization parameters and set up SR3 optimization for SINDy
threshold = 10^(-3);
kappa = 80; 
lambda = kappa*threshold^2/2;
stepsize = .01; 
trimmed_fraction = 0;
h = round((1-trimmed_fraction)*m);
maxiter = 100;
iter = 0;
err = 10;
dX = dXs(:,1:3);
w_old = A\dX;
v_old = ones(m,1);
% -- Solve SINDy regression problem using SR3
while err > 10^(-4) && iter < maxiter
    iter = iter + 1;
    Theta_c = A.*repmat(v_old,1,n);
    H = Theta_c.'*A+ kappa *eye(n); 
    % -- update coefficients
    coeffs = H\(Theta_c.'*dX + kappa*eye(n)*w_old);
    % -- update w by taking l_0 prox 
    w_new = coeffs;
    w_new(abs(w_new)<=threshold) = 0;
    % -- update v with gradient step and project onto capped h-simplex
    R2 = 0.5*sum((dX - A*coeffs).^2,2);
    v_new = proj_csimplex(v_old - stepsize*R2, h);
    
    % -- compute error
    err = norm(w_new - w_old)*kappa + norm(v_new-v_old)/stepsize;
    w_old = w_new; 
    v_old = v_new;
end

% -- plot resulting coefficients
c = lines;
nvars = 3;
figure
for j = 1:nvars
subplot(2,nvars,j)
plot(dX(:,j))
hold on
plot(A*coeffs(:,j))
subplot(2,nvars,j+3)
bar(coeffs(:,j),'FaceColor', c(j+3,:))
end
temp = coeffs;

dX = dXs(:,4:6);
A = [poolData(Xs(:,1:3),3,1,0), poolData(Xs(:,4:6),3,3,0)];
[m,n] = size(A);

% -- Set regularization parameters and set up SR3 optimization for SINDy
threshold = 4*10^(-5);
kappa = 80; 
lambda = kappa*threshold^2/2;
stepsize = .01; 
trimmed_fraction = 0;
h = round((1-trimmed_fraction)*m);
maxiter = 200;
iter = 0;
err = 10;

w_old = A\dX;
v_old = ones(m,1);
% -- Solve SINDy regression problem using SR3
while err > 10^(-8) && iter < maxiter
    iter = iter + 1;
    Theta_c = A.*repmat(v_old,1,n);
    H = Theta_c.'*A+ kappa *eye(n); 
    % -- update coefficients
    coeffs = H\(Theta_c.'*dX + kappa*eye(n)*w_old);
    % -- update w by taking l_0 prox 
    w_new = coeffs;
    w_new(abs(w_new)<=threshold) = 0;
    % -- update v with gradient step and project onto capped h-simplex
    R2 = 0.5*sum((dX - A*coeffs).^2,2);
    v_new = proj_csimplex(v_old - stepsize*R2, h);
    
    % -- compute error
    err = norm(w_new - w_old)*kappa +norm(v_new-v_old)/stepsize
    w_old = w_new; 
    v_old = v_new;
end

figure
for j = 1:nvars
subplot(2,nvars,j)
plot(dX(:,j))
hold on
plot(A*w_new(:,j))
subplot(2,nvars,j+3)
bar(w_new(:,j),'FaceColor', c(j,:))
end

size(coeffs);
% -- remove extra constant row
coeffs = [coeffs(1:4,:);coeffs(6:end,:)];
% -- pad first half
temp = [temp; zeros(16,3)];
coeffs = [temp,coeffs];

figure
imagesc(abs(coeffs))
caxis([0 0.1])
c = gray;
c = flipud(c);
colormap(c);
axis off
%% Simulate discovered dynamics
init = Xs(1,:).'; tspan = [0 1000]; 
[t,y] = ode45('DifferentLibraries_rhs',tspan,init,[],coeffs,3,1,0,3,0);
c = lines;
figure
for j = 1:6
    subplot(6,1,j)
    hold on
    plot(times(start:end)-times(start),X(:,j),'LineWidth',1.8, 'color', [c(j,:), 0.4])
    plot(t(:),y(:,j),':', 'LineWidth',3,'Color',c(j,:));
    xticklabels('')
    yticklabels('')  
end

%% Plot trimmed regions

A = poolData(Xs,1,3,0);
kappa = 80; 
stepsize = .01; 
trimmed_fraction = 0.2;
threshold =1-trimmed_fraction;
maxiter = 200;
[Xi_full,Xi_sparse,trimmed_points] = ...
    SR3(A, dXs,threshold, kappa, stepsize, trimmed_fraction, maxiter);
figure
hold on
plot3(x1,x2,x3,'k','LineWidth',1.8)
plot3(x1(trimmed_points<threshold),x2(trimmed_points<threshold),x3(trimmed_points<threshold), 'r.')
figure
hold on
plot3(x4,x5,x6,'k','LineWidth',1.8)
plot3(x4(trimmed_points<threshold),x5(trimmed_points<threshold),x6(trimmed_points<threshold), 'r.')