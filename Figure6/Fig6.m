% Figure 6
% Instead of traditional SINDy on the full trajectory of the FHN
% oscillator, we will try trimming and applying SINDy to each of the 6 resulting
% segments

% -- supporting files
% - FHN_rhs.m
% - poolData.m
% - SR3.m
% - cbrewer for Red-Gray divergent colormap

clear all; close all; clc
%% Simulate FHN network
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
t= t1(start:end)-t1(start);

x1 = Sr(1,1)*Vr(:,1); x4 = Sv(1,1)*Vv(:,1);
x2 = Sr(2,2)*Vr(:,2); x5 = Sv(2,2)*Vv(:,2);
x3 = Sr(3,3)*Vr(:,3); x6 = Sv(3,3)*Vv(:,3);
%% Apply SINDy with trimming to the second voltage PC
X = [ x1 x2 x3 x4 x5 x6];
Xs = X(2:end-1,:);
ts = t(2:end-1);
dt = t(2) - t(1);
dXdt = (X(3:end,:)-X(1:end-2,:))./(2*dt);

A = poolData(Xs,1,3,0);
kappa = 80; 
stepsize = .01; 
trimmed_fraction = 0.2;
threshold =1-trimmed_fraction;
maxiter = 200;
[Xi_full,Xi_sparse,trimmed_points] = ...
    SR3(A, dXdt(:,6),threshold, kappa, stepsize, trimmed_fraction, maxiter);

[pks, locs] = findpeaks(-trimmed_points,'MinPeakHeight',-.8);
delta = 40;
extended_trim = trimmed_points;
for j = 1:numel(locs)
    extended_trim(locs(j)-delta:locs(j)+delta) = 0;
end
figure
hold on
plot(-extended_trim)
plot(locs, zeros(size(locs)), 'ko')
inds = extended_trim;
inds(extended_trim<threshold) = 0;
inds(extended_trim>=threshold) = 1;

figure
subplot(1,2,1)
hold on
plot3(x1,x2,x3,'k','LineWidth',1.2)
plot3(x1(inds==0),x2(inds==0),x3(inds==0), 'r.')
grid on
subplot(1,2,2)
hold on
plot3(x4,x5,x6,'k','LineWidth',1.2)
plot3(x4(inds==0),x5(inds==0),x6(inds==0), 'r.')
grid on

figure
for j = 1:6
subplot(6,1,j)
plot(ts,Xs(:,j),'k')
hold on 
plot(ts(inds==0),Xs(inds==0,j),'r.')
set(gca, 'FontSize', 18,'FontName', 'Cambria')
end

%% Break trajectory into 6 subsections
inds1 = (inds == 1);
starts1 = [1;find((inds1(2:end) -inds1(1:end-1))==1)];
ends1 = [find((inds1(2:end) -inds1(1:end-1))==-1);length(inds)]; 
grp1a = zeros(size(inds));
    for j = 1:3:numel(starts1)
        s = starts1(j); f = ends1(j);
        grp1a(s:f) = 1;
    end
grp1b = zeros(size(inds));
    for j = 2:3:numel(starts1)
        s = starts1(j); f = ends1(j);
        grp1b(s:f) = 1;
    end
grp1c = zeros(size(inds));
    for j = 3:3:numel(starts1)
        s = starts1(j); f = ends1(j);
        grp1c(s:f) = 1;
    end
inds2 = (inds ==0);
starts2 = [find((inds2(2:end) -inds2(1:end-1))==1)];
ends2 = [find((inds2(2:end) -inds2(1:end-1))==-1)];

grp2a = zeros(size(inds));
    for j = 1:3:numel(starts2)
        s = starts2(j); f = ends2(j);
        grp2a(s:f) = 1;
    end

    grp2b = zeros(size(inds));
    for j = 2:3:numel(starts2)
        s = starts2(j); f = ends2(j);
        grp2b(s:f) = 1;
    end

    grp2c = zeros(size(inds));
    for j = 3:3:numel(starts2)
        s = starts2(j); f = ends2(j);
        grp2c(s:f) = 1;
    end
Grps = [grp1a, grp1b, grp1c, grp2a, grp2b, grp2c];
colors = lines(6);

for j = 1:6
ii = Grps(:,j)>0;
c = colors(j,:);
figure(17) 
hold on
plot3(Xs(ii>0,4),Xs(ii>0,5),Xs(ii>0,6),'.','color',c)
set(gca, 'FontSize', 18,'FontName', 'Cambria')
axis([-7,6, -5, 6,-2 2])
figure(18)
hold on
plot3(Xs(ii>0,1),Xs(ii>0,2),Xs(ii>0,3),'.','color',c)
end
set(gca, 'FontSize', 18,'FontName', 'Cambria')
%% Apply SINDy to each segment
taus = 10^(-3)*ones(1,6);
kappa = 80; 
stepsize = .01; 
trimmed_fraction = 0.1;
Xi = zeros(6,6,7);
for j = 1:6
    threshold = taus(j);
    ii = Grps(:,j) >0;
    x_ = Xs(ii,:);
    dx_ = dXdt(ii,:);
    A = poolData(x_,6,1,0);
    for jj = 1:6
        if jj<=3
            threshold = 1*10^(-3);
        end
        var = jj;
        [Xi_full,Xi_sparse,trimmed_points] = SR3(A, dx_(:,var),threshold,...
                                                kappa, stepsize, trimmed_fraction, 100);
        Xi(j,jj, :) = Xi_sparse;
    end
end
%% Visualize resulting fits against the derivative
figure
for j = 1:6
    inds = Grps(:,j) >0;
    x_ = Xs(inds,:);
    dx_ = dXdt(inds,:);
    A = poolData(x_,6,1,0);
    for jj = 1:6
    subplot(6,1,jj)
    hold on
    plot(t(inds), dx_(:,jj),'k.')
    plot(t(inds), A*squeeze(Xi(j,jj,:)),'.','color',colors(j,:))
    xticklabels('')
    yticklabels('')
    end
end

%% Visualize model coefficients
cmap =cbrewer('div', 'RdGy', 255, 'PCHIP');
figure
for j = 1:3
subplot(3,2,2*j-1)
imagesc([squeeze(Xi(j,1:3,:)*100).',squeeze(Xi(j,4:6,:)).'])
xticklabels('')
yticklabels('')
colormap(cmap)
caxis([-5 5])
colorbar
end
for j = 4:6
subplot(3,2,2*(j-3))
imagesc([squeeze(Xi(j,1:3,:)*1000).',squeeze(Xi(j,4:6,:)).'])
xticklabels('')
yticklabels('')
colormap(cmap)
caxis([-20 20])
colorbar
end
