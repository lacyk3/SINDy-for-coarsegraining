%% Run Simulations with different ratio of oscillators
clear all; close all; clc
%% Set up constant paramters

rng(447)
n = 100; % number of total oscillators

% -- Kuramoto parameters
a =.6;% range of possible frequencies is a - eps to a + eps
eps = 0.2;
KK = 10;

% -- FHN parameters 
alpha =[-.1,1.1,-1,.01,.01]; % parameters for self coupling function
KF =.5;  % coupling strength

opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
%% Loop through different nK values. 
% Run 1000 simulations and save the average dimension

Results = zeros(1000,100);
Results1 = zeros(1000,100);
Results2 = zeros(1000,100);
parfor (nK=1:100, 12) % number of Kuramoto oscillators, number of workers
    collectdims = zeros(1000,1);
    temp2 =zeros(1000,1);
        for j = 1:1000
            A=rand(n,n);  
            A=(A>0.8).*A; % adjacency matrix
            omega=a*(rand(nK,1)+0.5); % frequency for each oscillator
            nF = n - nK;
            thetai= a*2*randn(nK,1); % initial condition Kuramoto
            Vi= [-.95 + 2*rand(nF,1); -.05 + .1*rand(nF,1)]; % initial FHN
            x = [thetai;Vi];

            % -- time domain
            tspan =linspace(0,800,25000);
            dt = tspan(2) -tspan(1);

            [t,y] = ode45('HetNet_rhs',tspan,x,opts,nK,omega,KK,nF,alpha,KF,A);
            start = 10000;
            X = [cos(y(start:end,1:nK)) y(start:end,nK+1:n)].';

            [U,S,V] = svd(X,'econ');

            temp = find(sqrt(cumsum(diag(S).^2))./sqrt(sum(diag(S).^2))>.95);
            dim = temp(1);
            collectdims(j) = dim;
            
            temp = find(sqrt(cumsum(diag(S).^2))./sqrt(sum(diag(S).^2))>.99);
            dim = temp(1);
            temp2(j) = dim;
            
            temp = find(sqrt(cumsum(diag(S).^2))./sqrt(sum(diag(S).^2))>.9);
            dim = temp(1);
            temp3(j) = dim;
        end
        
    Results(:,nK) = temp3;
    Results1(:,nK) = collectdims;
    Results2(:,nK) = temp2;
    nK
end

save('RatioDimension90.mat','Results')
save('RatioDimension95.mat','Results1')
save('RatioDimension99.mat','Results2')

