% -- Rayleigh Network
% Each node has a voltage, v, and a recovery variable, w, following
% vi = wi
% eps wi = wi - wi^3/3 - vi + sum_j=1^n Aij(wj-wi)(1 + (vj-vi)^2)

% --supporting files
% - Rayleigh_rhs.m
% - poolData.m
% - sparsifyDynamics.m
% - sparseGalerkin.m
% - SR3.m

clear all; close all; clc
%% Plot behavior of one oscillator
eps = .001;
Rayleigh = @(t,x) ([x(2);
                    1/eps*(x(2) - x(2)^3/3 - x(1))]);
tspan = linspace(0,5,2000);
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45(Rayleigh,tspan,[0;-1.72],options);
figure
hold on
plot(t,y(:,1),'k','LineWidth',2)
plot(t,y(:,2),'k:','LineWidth',2)
axis off

%% -- Network of Rayleigh Oscillators w/ HKB coupling

% -- set up parameters
rng(288)
n=100;
eps = unifrnd(.02, .05, n,1);
A = ones(n,n);
K = .5;
on = 1;
x0 = randn(2*n,1);
tspan = linspace(0,15,2^17);

% -- simulate
options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
dt = tspan(2) - tspan(1);
[t0,x0] = ode45('Rayleigh_rhs',tspan,x0,options,eps,n,K,A,on);
[t,x] = ode45('Rayleigh_rhs',tspan,x0(end,:).',options,eps,n,K,A,on);

% -- raster plot
figure
h = pcolor(x(:,1:n).')
set(h,'EdgeColor','none');
colormap('jet')
colorbar()
xticklabels('')
set(gca, 'FontSize',18,'FontName','Cambria')
%% SVD
[Uv,Sv,Vv] = svd(x(:,1:n).','econ');
figure % voltage
plot(log(diag(Sv(1:10,1:10))/cumsum(diag(Sv(1:10,1:10)))),'ro')
figure
plot3(Sv(1,1)*Vv(:,1),Sv(2,2)*Vv(:,2),Sv(3,3)*Vv(:,3),'k','LineWidth',2);
grid on
set(gca, 'FontSize',18,'FontName','Cambria')

[Ur,Sr,Vr] = svd(x(:,n+1:2*n).','econ');
figure %recovery
plot(diag(Sr(1:10,1:10))./cumsum(diag(Sr(1:10,1:10))),'ro', 'LineWidth',1.5)
set(gca, 'FontSize',18,'FontName','Cambria')
figure 
plot3(Sr(1,1)*Vr(:,1),Sr(2,2)*Vr(:,2),Sr(3,3)*Vr(:,3),'k','LineWidth',2);
set(gca, 'FontSize',18,'FontName','Cambria')
figure %combo
plot3(Sv(1,1)*Vv(:,1),Sv(2,2)*Vv(:,2),Sr(2,2)*Vr(:,2),'k','LineWidth',2);

%% Try SINDy w/o trimming
X1 = zeros(size(x,1),2);
X1(:,1) = (Uv(:,1).'*x(:,1:n).').';
X1(:,2) = (Ur(:,1).'*x(:,n+1:2*n).').';
X1dot = (X1(3:end,:)-X1(1:end-2,:))/(2*dt);
A1 = poolData(X1(2:end-1,:),2,3,0);

X2 = zeros(size(x,1),2);
X2(:,1) = (Uv(:,2).'*x(:,1:n).').';
X2(:,2) = (Ur(:,2).'*x(:,n+1:2*n).').';
X2dot = (X2(3:end,:)-X2(1:end-2,:))/(2*dt);
A2 = poolData(X2(2:end-1,:),2,3,0);

% -- STLSQ to get coefficients
lam = 10^(-4);
coeffs1 = sparsifyDynamics(A1,X1dot(:,1),lam,1); 
coeffs2 = sparsifyDynamics(A1,X1dot(:,2),lam,1);
coeffs3 = sparsifyDynamics(A2,X2dot(:,1),lam,1); 
coeffs4 = sparsifyDynamics(A2,X2dot(:,2),lam,1);

options = odeset('RelTol',1e-7, 'AbsTol',1e-7);
Xi = [coeffs1,coeffs2];
init = X1(1,:).';
tspan = linspace(0,10, 20000); 
[t1,y1] = ode45('sparseGalerkin',tspan,init,options,Xi,3,0);
% 
figure
for j = 1:2
subplot(2,2,j)
plot(X1dot(:,j))
hold on
plot(A1*Xi(:,j), '.')
end
Xi = [coeffs3,coeffs4];
init = X2(1,:).';
[t2,y2] = ode45('sparseGalerkin',tspan,init,options,Xi,3,0);
for j = 1:2
subplot(2,2,2+j)
plot(X2dot(:,j))
hold on
plot(A2*Xi(:,j), '.')
end

% -- visualize coefficients
Xi_vis = zeros(12,4);
Xi_vis(1:3,1) = coeffs1(1:3);
Xi_vis(1:3,2) = coeffs2(1:3);
Xi_vis(1,3) = coeffs3(1);
Xi_vis(4:end,3) = coeffs3(2:end);
Xi_vis(1,4) = coeffs4(1);
Xi_vis(4:end,4) = coeffs4(2:end);
figure
imagesc(abs(Xi_vis))
caxis([0 16])
c = gray;
c = flipud(c);
colormap(c);
%set(gca, 'FontSize',18,'FontName','Cambria')
axis off
%% Simulate resulting model
colors = lines(6);
figure
subplot(4,1,1)
hold on
plot(t(1:10:85000),X1(1:10:85000,1),'color', [colors(1,:), 0.4],'LineWidth',1.8)
plot(t1(1:10:end),y1(1:10:end,1),':','color', colors(1,:),'LineWidth',2.8)
xticklabels('')
yticklabels('')
subplot(4,1,2)
hold on
plot(t(1:10:85000),X1(1:10:85000,2),'color', [colors(3,:), 0.4],'LineWidth',1.8)
plot(t1(1:10:end),y1(1:10:end,2),':','color', colors(3,:),'LineWidth',2.8)
xticklabels('')
yticklabels('')
subplot(4,1,3)
hold on
plot(t(1:10:85000),X2(1:10:85000,1),'color', [colors(2,:), 0.4],'LineWidth',1.8)
plot(t2(1:10:end),y2(1:10:end,1),':','color', colors(2,:),'LineWidth',2.8)
xticklabels('')
yticklabels('')
subplot(4,1,4)
hold on
plot(t(1:10:85000),X2(1:10:85000,2),'color', [colors(5,:), 0.4],'LineWidth',1.8)
plot(t2(1:10:end),y2(1:10:end,2),':','color', colors(5,:),'LineWidth',2.8)
set(gca, 'FontSize',18,'FontName','Cambria')
yticklabels('')

%% Where does trimming ID rapid dynamics?
kappa = 80; 
tau = 10^(-3);%threshold %10^(-4)
ss= .01; %step size
tf = 0.2;%trimmed fraction
maxiter = 600;

[coeffs_x1, Xi1, trimmed_points1] = SR3(A1,X1dot(:,1),tau,kappa, ss, tf, maxiter); 
[coeffs_x2, Xi2, trimmed_points2]  = SR3(A1,X1dot(:,2),tau,kappa, ss, tf, maxiter);
[coeffs_x3, Xi3, trimmed_points3]  = SR3(A2,X2dot(:,1),tau,kappa, ss, tf, maxiter); 
[coeffs_x4, Xi4, trimmed_points4]  = SR3(A2,X2dot(:,2),tau,kappa, ss, tf, maxiter); 

f = figure
hold on
plot3(X1(:,1),X1(:,2),X2(:,2),'k', 'LineWidth',2)
plot3(X1(trimmed_points4<.1,1),X1(trimmed_points4<.1,2),X2(trimmed_points4<.1,2),'r.','MarkerSize',15)
grid on
set(gca, 'FontSize',18,'FontName','Cambria')
view([-73.0667   17.4667])
%f.Position = [0.1300    0.1100    0.7750    0.8150]