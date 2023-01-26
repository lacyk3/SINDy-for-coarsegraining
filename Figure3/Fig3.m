%Fig 3
% SINDy with trimming to identify temporal boundary layers in limit cycle
% of Van der Pol oscillator. 
% The system of interest is: eps*y'' - (1-y'^2/3)y' + y = 0. 
       
% Making the substitution z = y' we get a system of two coupled ODEs:
% y' = z
% eps*z' = (1-y^2)z - y

% supporting function files:
% - SR3.m
% - Library.m (edit to match the BuildA() inline function)
% - hybrid_rhs.m
% - poolData.m
clear all;close all; clc
%% -- generate data
eps = 0.01;
VanderPol = @(t,x) [x(2);
                    1/eps*((1-x(1)^2)*x(2)- x(1))];
                
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45(VanderPol,[0 20],[-.8;-.1],options);
dx_centerdiff = (y(3:end,:) - y(1:end-2,:))./...
    [(t(3:end) - t(1:end-2)) (t(3:end) - t(1:end-2))];

%% Plot data, outer solution from approximation theory, and trimmed points
figure
hold on
plot([0,0],[-200,200],'k'); 
plot([-3,3],[0,0],'k')
plot(y(:,1),y(:,2),'k', 'LineWidth',2);
set(gca, 'FontSize', 18,'FontName', 'Cambria')
axis([-3 3 -200 200])

z_0 = linspace(-200,0,2000);
y_out_lower = (-1 - sqrt(1 + 4*z_0.^2))./(2*z_0);
plot(y_out_lower,z_0, 'k:','LineWidth',2); hold on

z_0 = linspace(0, 200,2000);
y_out_lower = (-1 - sqrt(1 + 4*z_0.^2))./(2*z_0);
plot(y_out_lower,z_0, 'k:','LineWidth',2); 

% -- SINDy via SR3 with trimming'
nvars = 2; polyorder = 3; usesine = 0;
A = poolData(y(2:end-1,:),nvars,polyorder,usesine);
tau = 1;
kappa = 80; 
% Note that SR3 calculates lambda = kappa*tau^2/2;
stepsize = .01; 
trimmed_fraction = 0.5;
threshold = 1-trimmed_fraction+.05; % boost 
[Xi_full,Xi_sparse,trimmed_points] = ...
    SR3(A, dx_centerdiff,tau, kappa, stepsize, trimmed_fraction, 2000);

% -- Visualize Results
plot(y(trimmed_points<threshold,1),y(trimmed_points<threshold,2),'r.', 'MarkerSize',10);
hold off

%% Construct partition from trimming results
ts = t(2:end-1);
ys = y(2:end-1,:);

% binarize trimmed_points
v  = trimmed_points;
v(trimmed_points < threshold) = 0;
v(trimmed_points >= threshold) = 1;

% difference v
v_diff = v(2:end) - v(1:end-1);
indices = 1:numel(ts);
% define start (v = -1) and end (v = 1) of trimmed intervals 
interval_bounds = [indices(v_diff == -1); indices(v_diff == 1)];

% Loop through trimmed intervals and make list of rectangular bounding
% boxes
rectangles = [];
figure
hold on
plot(ys(:,1),ys(:,2),'k.')
for j=1:size(interval_bounds,2)
  plot(y(interval_bounds(1,j):interval_bounds(2,j),1),...
      y(interval_bounds(1,j):interval_bounds(2,j),2),'b.') 
  plot(y(interval_bounds(1,j),1),y(interval_bounds(1,j),2),'g.')  
  plot(y(interval_bounds(2,j),1),y(interval_bounds(2,j),2),'r.') 
  x1 = min(y(interval_bounds(1,j):interval_bounds(2,j),1));
  x2 = max(y(interval_bounds(1,j):interval_bounds(2,j),1));
  y1 = min(y(interval_bounds(1,j):interval_bounds(2,j),2));
  y2 = max(y(interval_bounds(1,j):interval_bounds(2,j),2));
 
  h = rectangle('Position', [x1, y1, (x2-x1), y2-y1], ...
                'Curvature', 0.2, ...
                'FaceColor', [1, 0, 0, 0.1], ...
                'EdgeColor', [1, 0, 0, 0.1]);
            
  rectangles = [rectangles; [x1, y1, (x2-x1), y2-y1]];
end

partition = rectangles(1,:);

for j = 2:size(rectangles,1)
  % identify corners of new rectangle
  x1 = rectangles(j,1); 
  y1 = rectangles(j,2);
  x2 = rectangles(j,1)+rectangles(j,3);
  y2 = rectangles(j,2)+rectangles(j,4);
  corners = [x1, y1; x1, y2; x2,y2;x2,y1];
  % loop through rectangles already in the partition
  flag = 0;
  for p = 1:size(partition,1)     
      px1 = partition(p,1); 
      py1 = partition(p,2);
      px2 = partition(p,1)+ partition(p,3);
      py2 = partition(p,2)+ partition(p,4);
      pcorners = [px1, py1; 
                  px1, py2; 
                  px2,py2;
                  px2,py1];
      for c = 1:4
          % check whether corners lie in partition rectangle (or vice
          % versa)
          if (corners(c,1) >= px1) && (corners(c,1) <= px2) && ...
                  (corners(c,2) >= py1) && (corners(c,2) <= py2)
          % merge if so
          new_x1 = min(x1,px1); new_x2 = max(x2,px2);
          new_y1 = min(y1,py1); new_y2 = max(y2,py2);
          partition(p,:) = [new_x1, new_y1, new_x2-new_x1, new_y2-new_y1];
          flag = 1;
          elseif (pcorners(c,1) >= x1) && (pcorners(c,1) <= x2) && ...
                  (pcorners(c,2) >= y1) && (pcorners(c,2) <= y2)
          % merge if so
          new_x1 = min(x1,px1); new_x2 = max(x2,px2);
          new_y1 = min(y1,py1); new_y2 = max(y2,py2);
          partition(p,:) = [new_x1, new_y1, new_x2-new_x1, new_y2-new_y1];
          flag = 1;
          end
      end
  end
  if flag == 0
      % if not merged, append new rectangle to the existing partition.
      partition = [partition; rectangles(j,:)];
  end
end

figure
hold on
plot(ys(:,1),ys(:,2),'k.')
for i= 1:size(partition,1)
  h = rectangle('Position', partition(i,:), ...
                'Curvature', 0.2, ...
                'FaceColor', [1, 0, 0, 0.1], ...
                'EdgeColor', [1, 0, 0, 0.1]);
end

% Merge all in lower half plane and all in upper half plane
up = partition(partition(:,2)>0,:);
pup = [min(up(:,1)) min(up(:,2)) max(up(:,1)+up(:,3))-min(up(:,1)) max(up(:,2)+up(:,4))-min(up(:,2))];
down = partition(partition(:,2)<0,:);
pdown = [min(down(:,1)) min(down(:,2)) max(down(:,1)+down(:,3))-min(down(:,1)) max(down(:,2)+down(:,4))-min(down(:,2))];

figure
hold on
plot(ys(:,1),ys(:,2),'k.')
h = rectangle('Position', pup, ...
            'FaceColor', [1, 0, 0, 0.1], ...
            'EdgeColor', [1, 0, 0, 0.1]);
h = rectangle('Position', pdown, ...
            'FaceColor', [0, 0, 0, 0.1], ...
            'EdgeColor', [0, 0, 0, 0.1]);
        
partition = [pup;pdown];
partition(1,2) = -2/3;
partition(2,4) = 142.4681;
%% -- Translate partition on phase space to filters on input data

slowpoints = ones(length(ys),1);
fastpoints = zeros(length(ys),size(partition,1));

for i = 1:size(ys,1)
    x = ys(i,1); y = ys(i,2);
    for p = 1:size(partition,1)     
      px1 = partition(p,1); 
      py1 = partition(p,2);
      px2 = partition(p,1)+ partition(p,3);
      py2 = partition(p,2)+ partition(p,4);
      if (x >= px1) && (x <= px2) && ...
                  (y >= py1) && (y <= py2)
      slowpoints(i) = 0;
      fastpoints(i,p) = 1;
      end
    end
end

%% SINDy on dynamic regions

% -- Xslow
tslow = t(slowpoints>0);
Xslow =  ys(slowpoints>0,:);
dXslow = dx_centerdiff(slowpoints>0,:);
figure
hold on
plot(tslow,Xslow)
plot(tslow,dXslow)

X = Xslow(:,1); Y = Xslow(:,2);

% use SINDy with trimming to discount transition region
nvars = 1; polyorder = 3; usesine = 0;
A = [poolData(X,nvars,polyorder,usesine), X./(1-X.^2), X.*(1+X.^2)./(1-X.^2).^3];
labels = {'1','x','x^2','x^3', 'f1', 'f2'};
n = length(labels);

% hyperparameters
trimmed_fraction = .1; %0.2;
threshold = 1-trimmed_fraction;
tau = .8;
kappa = 80; 
stepsize = 1; 
maxiter = 2000;

% apply SR3 algorithm
[xXi_full,xXi_sparse,xtrimmed_points] = ...
    SR3(A, dXslow(:,1),tau, kappa, stepsize, trimmed_fraction, maxiter);

trimmed_fraction = .1; %0.2;
threshold = 1-trimmed_fraction;
tau = .8;
kappa = 80; 
stepsize = 1; 
maxiter = 2000;
[yXi_full,yXi_sparse,ytrimmed_points] = ...
    SR3(A, dXslow(:,2),tau, kappa, stepsize, trimmed_fraction, maxiter);
% -- plot
figure
subplot(2,2,1)
hold on
barh(xXi_sparse)
yticks(1:n)
yticklabels(labels)
subplot(2,2,2)
hold on
plot(tslow,dXslow(:,1),'.')
plot(tslow(xtrimmed_points<threshold),dXslow(xtrimmed_points<threshold,1),'.')
plot(tslow, A*xXi_sparse)
axis([1 5 -5 5])
subplot(2,2,3)
hold on
barh(yXi_sparse)
yticks(1:n)
yticklabels(labels)
subplot(2,2,4)
hold on
plot(tslow,dXslow(:,2),'.')
plot(tslow(ytrimmed_points<threshold),dXslow(ytrimmed_points<threshold,2),'.')
plot(tslow, A*yXi_sparse)
axis([1 5 -100 100])


%% -- Xfast
% crop off transient early behavior
% -- Xslow
tfast = t((fastpoints(:,1)+fastpoints(:,2))>0);
tfast = tfast(500:end);
Xfast =  ys((fastpoints(:,1)+fastpoints(:,2))>0,:);
Xfast = Xfast(500:end,:);
dXfast = dx_centerdiff((fastpoints(:,1)+fastpoints(:,2))>0,:);
dXfast = dXfast(500:end,:);
X = Xfast(:,1); Y = Xfast(:,2);
figure
subplot(2,1,1)
hold on
plot(tfast,dXfast(:,1),'.')
plot(tfast, Xfast(:,2),'.')
xlim([0 5])
subplot(2,1,2)
hold on
plot(tfast,dXfast(:,2),'.')
%plot(tfast, ,'.')
xlim([0 5])
%plot(tfast,dXfast,':')
%%
% SINDy w/o trimming

nvars = 1; polyorder = 5; usesine = 0;
A = [Y.^0 Y.*poolData(X,nvars,polyorder,usesine)];
labels = {'1','y','yx','yx^2','yx^3', 'yx^4','yx^5',};
n = length(labels);

lambda = 1;
iter = 0;
maxiter = 3000;
% compute Sparse regression via sequential least squares thresholding
wfasty = A\dXfast(:,2);
wfastx = A\dXfast(:,1);

while iter < maxiter
    smallinds = (abs(wfasty)<lambda);   % find small coefficients
    wfasty(smallinds)=0;                % and threshold           
    biginds = ~smallinds;
    wfasty(biginds) = A(:,biginds)\dXfast(:,2); 
    smallinds = (abs(wfastx)<.5);   % find small coefficients
    wfastx(smallinds)=0;                % and threshold           
    biginds = ~smallinds;
    wfastx(biginds) = A(:,biginds)\dXfast(:,1); 
    iter = iter + 1;
end

figure
subplot(2,2,1)
barh(wfasty)
yticks(1:n)
yticklabels(labels)
subplot(2,2,2)
hold on
plot(tfast,dXfast(:,2),'.')
plot(tfast, A*wfasty,'.')
xlim([1 5])
subplot(2,2,3)
barh(wfastx)
yticks(1:n)
yticklabels(labels)
subplot(2,2,4)
hold on
plot(tfast,dXfast(:,1),'.')
plot(tfast, A*wfastx,'.')
xlim([1 5])


%% Simulate results
% merge libraries
weights = [[zeros(6,1);wfastx]  [zeros(6,1);wfasty] ...
          [zeros(6,1);wfastx] [zeros(6,1);wfasty] ...
          [xXi_sparse; zeros(7,1)] [yXi_sparse; zeros(7,1)]]; 
init = Xfast(850,:).'; % start on fast segment
tspan = linspace(0,60,100000); 
%options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45('hybrid_rhs',tspan,init,[],weights,partition);
%%
figure
hold on 
plot(ys(:,1),ys(:,2))
plot(y(:,1),y(:,2),':')
plot(init(1),init(2),'o')

% Slow dynamics get around one corner, but not the other. Getting stuck in
% alternating between slow and fast around (-2.0142, -.0099). Our hybrid
% system has a little limit cycle there. Slow segment is also much too slow