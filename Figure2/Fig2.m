% Fig 2

% SINDy with trimming to identify temporal boundary layers in limit cycle
% of Rayleigh oscillator. The system of interest is: eps*y'' - (1-y'^2/3)y' + y = 0. 
       
% Making the substitution z = y' we get a system of two coupled ODEs:
% y' = z
% eps*z' = -y + (1-z^2/3)z

% supporting function files:
% - Library.m (edit to match the BuildA() inline function)
% - hybrid_rhs.m

clear all; close all; clc
%% Simulate Rayleigh oscillator to create "data"

eps = .001;
init = [-.6771;-1];
Rayleigh = @(t,x) ([x(2);
                    1/eps*(x(2) - x(2)^3/3 - x(1))]);
                
tspan = linspace(0,2.5,10000);
options = odeset('RelTol',1e-8, 'AbsTol',1e-8);
[t,y] = ode45(Rayleigh,tspan,init,options);

dx_centerdiff = (y(3:end,:) - y(1:end-2,:))./...
                [(t(3:end) - t(1:end-2)) (t(3:end) - t(1:end-2))];
%% (1) SINDy with trimming
% build library 
nvars = 1; polyorder = 5; usesine = 0;
A = poolData(y(2:end-1,2),nvars,polyorder,usesine);
% hyperparameters
trimmed_fraction = .011; %0.2;
lambda = .989; 
kappa = 80; 
stepsize = .01; 
maxiter = 2000;

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
%% (2) Partition phase space
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

% Loop through list of rectangles and form partition by merging overlapping. 
% Note that this procedure only works because our regions of fast dynamics
% that we might want to treat separately are far apart, so there will only
% ever be one rectangle in the partition that we want to merge with at a
% time. Furthermore, growing the rectangles by merging them will not make 2
% rectangles that are already in the partition overlap, so we do not need to 
% re-check the partition This algorithm could be made more robust at a later date.

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
          % check whether corners lie in partition rectangle 
          if (corners(c,1) >= px1) && (corners(c,1) <= px2) && ...
                  (corners(c,2) >= py1) && (corners(c,2) <= py2)
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

%% (3) SINDy to discover dynamics on each section

BuildA = @(X,Y)[ X.^0, X, Y,  X.^2, Y.^2,...
                Y./((1./(Y-1).^2)+(2/3)*(1./(Y-1))-(2/3)*(1./(Y+2))),...
                Y./((1./(Y+1).^2)-(2/3)*(1./(Y+1))+(2/3)*(1./(Y-2))),...
                log(abs(X-1)), log(abs(X+1)),  ...
                1./(X-1), 1./(X+1), X./(1-X.^2),...
                Y./(1-Y.^2)];

% -- X_slow
Xslow =  ys(trimmed_points>=threshold,:);
dXslow = dx_centerdiff(trimmed_points>=threshold,:);

lambda = 10^(-1);
iter = 0;
maxiter = 20;

Aslowx = BuildA(Xslow(:,1),Xslow(:,2));
% compute Sparse regression via sequential least squares thresholding
wslowx = Aslowx\dXslow(:,1);  % initial guess: Least-squares
while iter < maxiter
    smallinds = (abs(wslowx)<lambda);   % find small coefficients
    wslowx(smallinds)=0;                % and threshold           
    biginds = ~smallinds;
    % Regress dynamics onto remaining terms to find sparse Xi
    wslowx(biginds) = Aslowx(:,biginds)\dXslow(:,1);   
    iter = iter + 1;
end

figure
barh(wslowx)
yticks(1:13)
yticklabels({'1','x','y','x^2','y^2','f1(y)','f2(y)','log|x-1|','log|x+1|',...
            '1/(x-1)','1/(x+1)', 'x/(1-x.^2)',...
             'y/(1-y^2)'})

% -- filter down to center of slow intervals         
Aslowy = BuildA([Xslow(1000:2000,1); ...
            Xslow(4000:5000,1);...
            Xslow(7000:8000,1)],...
            [Xslow(1000:2000,2); ...
            Xslow(4000:5000,2);...
            Xslow(7000:8000,2)]);
lambda = .8;
iter = 0;
maxiter = 10;

% compute Sparse regression via sequential least squares thresholding
wslowy = Aslowy\[dXslow(1000:2000,2); ...
            dXslow(4000:5000,2);...
            dXslow(7000:8000,2)];  % initial guess: Least-squares
wslowy(1) = 0; wslowy(4) = 0; wslowy(5) = 0; % and zero out unwanted terms
while iter < maxiter
    smallinds = (abs(wslowy)<lambda);   % find small coefficients
    wslowy(smallinds)=0;                % and threshold           
    biginds = ~smallinds;
    % Regress dynamics onto remaining terms to find sparse Xi
    wslowy(biginds) = Aslowy(:,biginds)\[dXslow(1000:2000,2); ...
            dXslow(4000:5000,2);...
            dXslow(7000:8000,2)];   
    iter = iter + 1;
end

figure
barh(wslowy)
yticks(1:13)
yticklabels({'1','x','y','x^2','y^2', 'f1(y)','f2(y)','log|x-1|','log|x+1|',...
           '1/(x-1)','1/(x+1)', 'x/(1-x.^2)',...
            'y/(1-y^2)'})

        
%% -- X_fast
Xfast =  ys(trimmed_points<threshold,:);
dXfast = dx_centerdiff(trimmed_points<threshold,:);
tfast = ts(trimmed_points<threshold);
lambda = 100;
iter = 0;
maxiter = 2000;

Afasty = BuildA(Xfast(:,1),Xfast(:,2));
% compute Sparse regression via sequential least squares thresholding
wfasty = Afasty\dXfast(:,2);
wfasty(8) = 0; wfasty(9) = 0; wfasty(10) = 0; wfasty(11) = 0;
while iter < maxiter
    smallinds = (abs(wfasty)<lambda);   % find small coefficients
    wfasty(smallinds)=0;                % and threshold           
    biginds = ~smallinds;
    % Regress dynamics onto remaining terms to find sparse Xi
    %wfasty(biginds) = Afasty(:,biginds)\[dXfast(1:36,2);dXfast(75:100,2)];  
    wfasty(biginds) = Afasty(:,biginds)\dXfast(:,2); 
    iter = iter + 1;
end

figure
barh(wfasty)
yticks(1:13)
yticklabels({'1','x','y','x^2','y^2','f1(y)','f2(y)','log|x-1|','log|x+1|',...
            '1/(x-1)','1/(x+1)', 'x/(1-x.^2)',...
             'y/(1-y^2)'})
figure
hold on
plot(tfast,Afasty*wfasty);
plot(tfast,dXfast(:,2),'.')

f1 = @(x) -(1/eps)*x./((1./(x-1).^2)+(2/3)*(1./(x-1))-(2/3)*(1./(x+2)));
f2 = @(x) -(1/eps)*x./((1./(x+1).^2)-(2/3)*(1./(x+1))+(2/3)*(1./(x-2)));


%% (4) Simulate system
% concatenate two sets of fast dynamics (one for each boundary layer) + slow dynamics 
weights = [zeros(13,1) wfasty zeros(13,1) wfasty wslowx wslowy]; 
init = [0;sqrt(3)]; % start on slow segment
tspan = linspace(0,2,1000); 
[t,y] = ode45('hybrid_rhs',tspan,init,[],weights,partition);
figure
plot(t,y)

% Mess with time step/ IC to make simulation stable for full time
