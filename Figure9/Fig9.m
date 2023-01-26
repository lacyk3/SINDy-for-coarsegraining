% -- Figure 9
% Simulate, coarse-grain, and trim networks of coupled
% Kuramoto, FHN and Rossler oscillators.

% -- supporting files
% - HetNet3.m
% - SR3.m
% - poolData.m
% - sparseGalerkin.m
clear all; close all; clc;
%% Set network wide parameters
rng(72)
n = 100; % number of total oscillators
A=rand(n,n);  % connectivity matrix
A=(A>0.7).*A; 
% -- Kuramoto nodes
a =.4;% range of possible frequencies is 0 to a, centered at a/2
eps = 0.2;
KK = 20; % coupling strength
% -- FHN nodes
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5; % coupling strength
% -- Rossler nodes
params = [0.2, 0.2, 5.7]; % parameters for self coupling function
KR = 10;  % coupling strength

%%  Row (a): Big, synchronized spiking behavior
% nK= 85; % number of Kuramoto oscillators
% omega= a + eps*(rand(nK,1)-0.5); % random natural frequency for each Kura oscillator
% nF = 10; 
% nR = n - nK - nF;
% %-- initial conditions
% thetai= a*2*randn(nK,1); % initial condition
% Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initiate state
% Ri = 2*randn(3*nR,1);
% x = [thetai;Vi;Ri];
% %-- simulate
% tspan =linspace(0,2000,27000);
% dt = tspan(2) -tspan(1);
% opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
% [t,y] = ode45('HetNet3',tspan,x,opts,nK,omega,KK,nF,alpha,KF,nR,params,KR,A, 1);
% %-- raster plot 
% f = figure
% h = pcolor([cos(y(:,1:nK)) y(:,nK+1:nK+2*nF) y(:,nK+2*nF+1:end)/10].');
% xticklabels('')
% yticklabels('')
% set(h, 'EdgeColor', 'none');
% colormap('jet')
% caxis([-1 1]);
% colorbar
% set(gca, 'FontSize', 18,'FontName', 'Cambria')
% f.Position = [100 100 400 300]
% %-- zoom in on last quarter of time
% start = 20000;
% f = figure
% h = pcolor([cos(y(start:end,1:nK)) y(start:end,nK+1:nK+2*nF) y(start:end,nK+2*nF+1:end)/10].');
% xticklabels('')
% yticklabels('')
% set(h, 'EdgeColor', 'none');
% colormap('jet')
% caxis([-1 1]);
% colorbar
% set(gca, 'FontSize', 18,'FontName', 'Cambria')
% f.Position = [100 100 400 300]
% 
% %-- low dimensional dynamics
% X = [cos(y(:,1:nK)) y(:,nK+1:nK+2*nF) y(:,nK+2*nF+1:end)/10].';
% [U,S,V] = svd(X(:,start:end),'econ');
% f = figure
% hold on
% plot3(U(:,3).'*X(:,start:end),U(:,4).'*X(:,start:end), U(:,1).'*X(:,start:end),'k', 'LineWidth',2)
% grid on
% view([ -211.2000   33.4667])
% set(gca, 'FontSize',18, 'FontName','Cambria')
% f.Position = [100 100 400 300]
% 
% %-- simulate at higher resolution from end of previous 
% tspan =linspace(0,270,2^17);
% dt = tspan(2) -tspan(1);
% opts = odeset('RelTol',1e-8, 'AbsTol',1e-8);
% [T,Y] = ode45('HetNet3',tspan,y(end,:).',opts,nK,omega,KK,nF,alpha,KF,nR,params,KR,A, 1);
% X = [cos(Y(:,1:nK)) Y(:,nK+1:nK+2*nF) Y(:,nK+2*nF+1:end)/10].';
% 
% 
% %%-- trim
% Xs = [(U(:,1).'*X).', (U(:,2).'*X).',(U(:,3).'*X).',(U(:,4).'*X).'];
% dxdt= (Xs(3:end,2)-Xs(1:end-2,2))/(2*dt);
% nvars = 4; polyorder = 3; usesine = 0;
% A = poolData(Xs(2:end-1,:),nvars,polyorder,usesine);
% kappa = 80; 
% tau = 10^(-5);%threshold 
% ss= .01; %step size
% tf = 0.2;%trimmed fraction
% maxiter = 600;
% [Xi_full,Xi_sparse,trimmed_points] = SR3(A, dxdt,tau, kappa, ss, tf, maxiter);
% 
% % -- These dynamics are noisier, so trimming doesn't just pull out regions
% % of rapid change
% threshold = 1-tf-;
% figure
% hold on
% plot3(Xtemp(:,3),Xtemp(:,4), Xtemp(:,1),'k.');
% plot3(Xtemp(find(trimmed_points<threshold),3),...
%     Xtemp(find(trimmed_points<threshold),4), ...
%     Xtemp(find(trimmed_points<threshold),1),'rx');
% view([ -211.2000   33.4667])
% 
% % -- Let's pull out the peaks from trimming
% [pks, locs] = findpeaks(-trimmed_points,'MinPeakHeight',-.74,'MinPeakDistance',500);
% figure
% hold on
% plot(-trimmed_points)
% plot(locs,pks, 'o')
% % and put a window around those peaks
% delta = 450;
% plot([locs-delta locs+delta],[pks pks] , 'k+')
% extended_trim = zeros(size(trimmed_points));
% for j = 1:numel(locs)
%     extended_trim(locs(j)-delta:locs(j)+delta) = 1;
% end
% % This looks like it is narrowing our focus to regions of rapid change
% f = figure
% hold on
% plot3(Xtemp(:,3),Xtemp(:,4), Xtemp(:,1),'k.');
% plot3(Xtemp(extended_trim>0,3),Xtemp(extended_trim>0,4), Xtemp(extended_trim>0,1),'rx');
% view([ -211.2000   33.4667])
% set(gca, 'FontSize',18, 'FontName','Cambria')
% grid on
% f.Position = [100 100 400 300]
% 
% f = figure
% for j = 1:4
% subplot(4,1,j)
% hold on
% plot(th,Xtemp(:,j),'k', 'LineWidth', 2)
% plot(th(extended_trim>0),Xtemp(extended_trim>0,j),'rx')
% xticklabels('')
% yticklabels('')
% box off
% end
% f.Position = [100 100 400 300]
%% Row (b): FHN dominant
% nK= 20; 
% omega= a + eps*(rand(nK,1)-0.5); 
% nF = 78; 
% nR = n - nK - nF;
% %-- initial conditions
% thetai= a*2*randn(nK,1); 
% Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; 
% Ri = 2*randn(3*nR,1);
% x = [thetai;Vi;Ri];
% %-- simulate
% tspan =linspace(0,2100,27000);
% dt = tspan(2) -tspan(1);
% opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
% [t,y] = ode45('HetNet3',tspan,x,opts,nK,omega,KK,nF,alpha,KF,nR,params,KR,A, 1);
% 
% %-- raster plot 
% f = figure
% h = pcolor([cos(y(:,1:nK)) y(:,nK+1:nK+2*nF) y(:,nK+2*nF+1:end)/10].');
% xticklabels('')
% yticklabels('')
% set(h, 'EdgeColor', 'none');
% colormap('jet')
% caxis([-1 1]);
% colorbar
% set(gca, 'FontSize', 18,'FontName', 'Cambria')
% f.Position = [100 100 400 300]
% 
% %-- low dimensional dynamics
% X = [cos(y(:,1:nK)) y(:,nK+1:nK+2*nF) y(:,nK+2*nF+1:end)/10].';
% [U,S,V] = svd(X(:,start:end),'econ');
% 
% f = figure
% hold on
% plot3(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end), U(:,3).'*X(:,start:end),'k', 'LineWidth',2)
% grid on
% view([ 125.6000   14.8000])
% set(gca, 'FontSize',18, 'FontName','Cambria')
% f.Position = [100 100 400 300]
% 
% %-- time series 
% % Use last state as initial condition and simulate at high resolution for a
% % short time
% Init = y(end,:).';
% tspan =linspace(0,400,2^17);
% dt = tspan(2) -tspan(1);
% opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
% [t2,y2] = ode45('HetNet3',tspan,Init,opts,nK,omega,KK,nF,alpha,KF,nR,params,KR,A, 1);
% 
% f = figure
% h = pcolor([cos(y2(:,1:nK)) y2(:,nK+1:nK+2*nF) y2(:,nK+2*nF+1:end)/10].');
% xticklabels('')
% yticklabels('')
% set(h, 'EdgeColor', 'none');
% colormap('jet')
% caxis([-1 1]);
% colorbar
% set(gca, 'FontSize', 18,'FontName', 'Cambria')
% f.Position = [100 100 400 300]
% X2 = [cos(y2(:,1:nK)) y2(:,nK+1:nK+2*nF) y2(:,nK+2*nF+1:end)/10].';

%%-- trim on lowrank dynamics
% [U,S,V] = svd(X2,'econ');
% Xs = [(U(:,1).'*X2).', (U(:,2).'*X2).',(U(:,3).'*X2).',(U(:,4).'*X2).'];
% dxdt= (Xs(3:end,1)-Xs(1:end-2,1))/(2*dt);
% nvars = 4; polyorder = 4; usesine = 0;
% A = poolData(Xs(2:end-1,:),nvars,polyorder,usesine);
% kappa = 80; 
% tau = 10^(-5);%threshold
% ss= .01; %step size
% tf = 0.2;%trimmed fraction
% maxiter = 800;
% [Xi_full,Xi_sparse,trimmed_points] = SR3(A, dxdt,tau, kappa, ss, tf, maxiter);
% 
% % straight trimming does okay, but sort of noisy
% figure
% m = 1 - tf;
% Xtemp = Xs(2:end-1,:);
% th = t2(2:end-1);
% subplot(2,1,1)
% plot(th.',Xs(2:end-1,:)); hold on
% plot(th(find(trimmed_points<m)).',Xtemp(find(trimmed_points<m),:),'rx')
% set(gca, 'FontSize',18)
% subplot(2,1,2)
% plot(th.',trimmed_points)
% 
% figure
% hold on
% plot3(Xtemp(:,2),Xtemp(:,3), Xtemp(:,4),'k');
% (plot3(Xtemp(find(trimmed_points<m),2), ...
%       Xtemp(find(trimmed_points<m),3), ...
%       Xtemp(find(trimmed_points<m),4),'rx'));
%   
% %%-- Try peak detection + window
% [pks, locs] = findpeaks(-trimmed_points,'MinPeakHeight',-.7995,'MinPeakDistance',5000);
% figure
% hold on
% plot(-trimmed_points)
% plot(locs,pks, 'o')
% delta = 1500;
% plot([locs-delta locs+delta],[pks pks] , 'k+')
% 
% extended_trim = zeros(size(trimmed_points));
% for j = 1:numel(locs)
%     extended_trim(max(locs(j)-delta,1):min(locs(j)+delta,numel(extended_trim)) )= 1;
% end
% % -- fill in the small gaps
% extended_trim(locs(14)-delta:locs(15)+delta) =1;
% 
% f=figure
% hold on
% plot3(Xtemp(:,2),Xtemp(:,3), Xtemp(:,1),'k.');
% plot3(Xtemp(extended_trim>0,2),Xtemp(extended_trim>0,3), Xtemp(extended_trim>0,1),'rx');
% view([ -211.2000   33.4667])
% f.Position = [100 100 400 300]
% set(gca, 'FontSize',18, 'FontName','Cambria')
% grid on
% % 
% f = figure
% for j = 1:4
% subplot(4,1,j)
% hold on
% plot(th,Xtemp(:,j),'k', 'LineWidth', 2)
% plot(th(extended_trim>0),Xtemp(extended_trim>0,j),'rx')
% xticklabels('')
% yticklabels('')
% box off
% end
% f.Position = [100 100 400 300]

%% row (c) Small oscillations
rng(2024) 
nK= 78; 
omega= a + eps*(rand(nK,1)-0.5); 
nF = 15; 
nR = n - nK - nF;
%-- initial conditions
thetai= a*2*randn(nK,1);
Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; 
Ri = 2*randn(3*nR,1);
x = [thetai;Vi;Ri];
%-- simulate
tspan =linspace(0,2000,27000);
dt = tspan(2) -tspan(1);
opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
[t,y] = ode45('HetNet3',tspan,x,opts,nK,omega,KK,nF,alpha,KF,nR,params,KR,A, 1);

%-- raster plot 
f = figure
h = pcolor([cos(y(:,1:nK)) y(:,nK+1:nK+2*nF) y(:,nK+2*nF+1:end)/10].');
xticklabels('')
yticklabels('')
set(h, 'EdgeColor', 'none');
colormap('jet')
caxis([-1 1]);
colorbar
set(gca, 'FontSize', 18,'FontName', 'Cambria')
f.Position=([100,100,400,300])

%-- low dimensional dynamics
start = 2000
X = [cos(y(:,1:nK)) y(:,nK+1:nK+2*nF) y(:,nK+2*nF+1:end)/10].';
[U,S,V] = svd(X(:,start:end),'econ');

f = figure
hold on
plot3(U(:,1).'*X(:,start:end),U(:,2).'*X(:,start:end), U(:,3).'*X(:,start:end),'k', 'LineWidth',2)
grid on
view([ 136.0000   19.6000])
set(gca, 'FontSize',18, 'FontName','Cambria')
f.Position=([100 100 400 300])

%-- trim
t1 = 20000;
Xs = [(U(:,1).'*X(:,t1:end -2600)).', (U(:,2).'*X(:,t1:end -2600)).',(U(:,3).'*X(:,t1:end -2600)).',(U(:,4).'*X(:,t1:end -2600)).'];
dxdt= (Xs(3:end,1)-Xs(1:end-2,1))/(2*dt);
nvars = 4; polyorder = 3; usesine = 0;
A = poolData(Xs(2:end-1,:),nvars,polyorder,usesine);
kappa = 80; 
tau = .1;%threshold %10^(-4)
ss= .01; %step size
tf = 0.2;%trimmed fraction
maxiter = 600;
[Xi_full,Xi_sparse,trimmed_points] = SR3(A, dxdt,tau, kappa, ss, tf, maxiter);

[pks, locs] = findpeaks(-trimmed_points,'MinPeakHeight',-.79995);
figure
hold on
plot(-trimmed_points)
plot(locs,pks, 'o')
delta = 6;
plot([locs-delta locs+delta],[pks pks] , 'k+')

extended_trim = zeros(size(trimmed_points));
for j = 1:numel(locs)
    extended_trim(locs(j)-delta:locs(j)+delta)= 1;
end

Xtemp = Xs(2:end-1,:);
th = t(t1:t1+length(Xtemp(:,1))-1)-t(t1);
f = figure
hold on
plot3(Xtemp(:,1),Xtemp(:,2), Xtemp(:,3),'k.');
plot3(Xtemp((extended_trim>0),1),Xtemp((extended_trim>0),2), Xtemp((extended_trim>0),3),'rx');
grid on
view([ 136.0000   19.6000])
set(gca, 'FontSize',18, 'FontName','Cambria')
f.Position=([100 100 400 300])

f = figure
for j = 1:3
subplot(3,1,j)
hold on
plot(th,Xtemp(:,j),'k', 'LineWidth', 2)
plot(th(extended_trim>0).',Xtemp(extended_trim>0,j),'rx');
xticklabels('')
yticklabels('')
box off
end
f.Position=([100 100 400 300])

%%--Try fitting dynamics and running result
Xs = [(U(:,1).'*X(:,t1:end)).',(U(:,2).'*X(:,t1:end)).',(U(:,3).'*X(:,t1:end)).'];

dxdt= (Xs(3:end,:)-Xs(1:end-2,:))/(2*dt);

nvars = 3; polyorder = 3; usesine = 1;
A = poolData(Xs(2:end-1,:),nvars,polyorder,usesine);
kappa = 80; 
tau = 10^(-5);
ss= .01; %step size
tf = 0.2;%trimmed fraction
maxiter = 600;
[Xi_full,Xi_sparse,trimmed_points] = SR3(A, dxdt,tau, kappa, ss, tf, maxiter);

opts = odeset('RelTol',1e-8, 'AbsTol',1e-8);
init = Xs(1,:).';
tspan = [0 100]; 
[ty,y] = ode45('sparseGalerkin',tspan,init,opts,Xi_sparse,polyorder,usesine);
figure
hold on
plot3(Xs(:,1),Xs(:,2),Xs(:,3))
plot3(y(:,1),y(:,2),y(:,3))