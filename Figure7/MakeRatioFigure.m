%% Make Figure Exploring Oscillator Ratio 
clear all; close all; clc
%% Changing Ratio of Kuramoto to FHN oscillators
f = figure
hold on
x = 0:99;
for j=[30,53,65,80]
    plot([j j], [-2 16],'k--', 'LineWidth',1.2)
end

% -- 99%
load('RatioDimension99.mat')
y1 = mean(Results2)+std(Results2);
y2 = mean(Results2)-std(Results2);
patch([x fliplr(x)], [y1 fliplr(y2)],[0.7 0.7 .7],'EdgeColor',[.7 0.7 .7])
alpha(0.3)

% -- 90%
load('RatioDimension90.mat')
y1 = mean(Results)+std(Results);
y2 = mean(Results)-std(Results);
patch([x fliplr(x)], [y1 fliplr(y2)],[0.4660 0.6740 0.1880],'EdgeColor',[0.4660 0.6740 0.1880])
alpha(0.25)


% -- 95%
load('RatioDimension95.mat')
y1 = mean(Results)+std(Results);
y2 = mean(Results)-std(Results);
patch([x fliplr(x)], [y1 fliplr(y2)],[0.3010 0.7450 0.9330],'EdgeColor',[0.3010 0.7450 0.9330])
alpha(0.3)



load('RatioDimension99.mat')
plot(0:99, mean(Results2), 'k','LineWidth',2)
load('RatioDimension95.mat')
plot(0:99, mean(Results),'color',[0.3010 0.7450 0.9330],'LineWidth',2)
load('RatioDimension90.mat')
plot(0:99, mean(Results), 'color',[0.4660 0.6740 0.1880],'LineWidth',2)
set(gca,'FontSize',18)
set(gca, 'FontName','Cambria')
% xticklabels('')
% yticklabels('')
box off
axis([0 99 -2 16])
f.Position=([10,100,1300,300])

%% Simulated Example at 30
rng(72) 
n = 100; % number of total oscillators
A=rand(n,n);  
A=(A>0.8).*A; % connectivity

% -- Kuramoto parameters
a =.6;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 10;
nK=30; % number of Kuramoto oscillators
omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each oscillator
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;
%-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];


% -- time domain
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);

opts = odeset('RelTol',1e-8,'AbsTol',1e-5);

[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);

f1 = figure
h = pcolor([cos(y(:,1:nK)) y(:,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f1.Position = [100   100    400   300]

f2 = figure
h = pcolor([cos(y(end-5000:end,1:nK)) y(end-5000:end,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f2.Position = [100   100    400    300]

start = 10000;

X = [cos(y(:,1:nK)) y(:,nK+1:n)].';

[U,S,V] = svd(X(:,start:end),'econ');

start = 13000;

f30=figure
plot3(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end), U(:,3).'*X(:,start:end),'k', 'LineWidth',1.5)
set(gca, 'FontSize',20, 'FontName', 'Cambria')
grid on
view([-105.3333   35.0667])
f30.Position = [100   100    300    300]
%% Simulated Example at 53
rng(1024)

n = 100; % number of total oscillators
A=rand(n,n);  
A=(A>0.8).*A; % connectivity

% -- Kuramoto parameters
a =.6;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 10;
nK=53; % number of Kuramoto oscillators
omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each oscillator
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;
%-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];


% -- time domain
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);

opts = odeset('RelTol',1e-8,'AbsTol',1e-5);

[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);

f1 = figure
h = pcolor([cos(y(:,1:nK)) y(:,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f1.Position = [100   100    400   300]

f2 = figure
h = pcolor([cos(y(end-5000:end,1:nK)) y(end-5000:end,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f2.Position = [100   100    400    300]

start = 10000;

X = [cos(y(:,1:nK)) y(:,nK+1:n)].';

[U,S,V] = svd(X(:,start:end),'econ');

start = 1;

f53 = figure;
hold on
plot3(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end), U(:,3).'*X(:,start:end),'k', 'LineWidth',1.5)
plot3(U(:,1).'*X(:,start),U(:,2).'*X(:,start), U(:,3).'*X(:,start),'.','color',[0.4660 0.6740 0.1880], 'MarkerSize',20)
plot3(U(:,1).'*X(:,end),U(:,2).'*X(:,end), U(:,3).'*X(:,end),'r.', 'MarkerSize',20)
set(gca, 'FontSize',20, 'FontName', 'Cambria')
grid on
f53.Position = [100    100    300    300];
view([-49.6000   27.0667]);

%% Simulated Example at 65
rng(30) 
%rng(200)

n = 100; % number of total oscillators
A=rand(n,n);  
A=(A>0.8).*A; % connectivity

% -- Kuramoto parameters
a =.6;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 10;
nK=65; % number of Kuramoto oscillators
omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each oscillator
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;
%-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];


% -- time domain
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);

opts = odeset('RelTol',1e-8,'AbsTol',1e-5);

[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);

f1 = figure
h = pcolor([cos(y(:,1:nK)) y(:,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f1.Position = [100   100    400   300]

f2 = figure
h = pcolor([cos(y(end-5000:end,1:nK)) y(end-5000:end,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f2.Position = [100   100    400    300]

start = 10000;

X = [cos(y(:,1:nK)) y(:,nK+1:n)].';

[U,S,V] = svd(X(:,start:end),'econ');

f65=figure
plot3(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end), U(:,3).'*X(:,start:end),'k', 'LineWidth',1.5)
set(gca, 'FontSize',20, 'FontName', 'Cambria')
grid on
view([-64.8000    4.6667])
%view([-22.6667   59.0667])
f65.Position = [100   100    300    300]
%% Simulated Example at 80
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
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength
nF = n - nK;
%-- initial conditions
thetai= a*2*randn(nK,1); % initial condition
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
x = [thetai;Vi];


% -- time domain
tspan =linspace(0,1000,20000);
dt = tspan(2) -tspan(1);

opts = odeset('RelTol',1e-8,'AbsTol',1e-5);

[t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);

f1 = figure
h = pcolor([cos(y(:,1:nK)) y(:,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f1.Position = [100   100    400   300]

f2 = figure
h = pcolor([cos(y(end-5000:end,1:nK)) y(end-5000:end,nK+1:n)].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f2.Position = [100   100    400    300]

start = 10000;

X = [cos(y(:,1:nK)) y(:,nK+1:n)].';

[U,S,V] = svd(X(:,start:end),'econ');

start = 2000
f80=figure
plot3(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end), U(:,3).'*X(:,start:end),'k', 'LineWidth',1.5)
set(gca, 'FontSize',20, 'FontName', 'Cambria')
grid on
view([ -39.4667   15.8667])
f80.Position = [100   100    300    300]