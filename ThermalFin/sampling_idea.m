%%
% The basic idea is to see if we can get an idea of where to sample the
% gradient by looking at the singular vectors of the gradient samples
% corresponding to spatial points. How will we look at them? We'll project
% their inputs using the same regression graphics approach. 
%
% Let's see what happens!

clear all
% Tan Bui, Mar 03 2014
% Files for Paul

Globals1D

al = 0.6;
alpha_a = 4;
alpha_b = -1;
alpha = alpha_a * 10^alpha_b;
elements = 32;
noise_a = 1;
noise_b = -2;
noise = noise_a * 10^noise_b;

%%
% This file contains the MAP point and the mesh structure, etc
% GamMap is the MAP point
st = strcat('PoissonMAP1D',num2str(elements),'al',num2str(al),'alpha',num2str(alpha),'noise',num2str(noise),'.mat');
load(st,'GamMap', 'data', 'mesh');

GeometryOrder = mesh.GeometryOrder; 
SolutionOrder = mesh.SolutionOrder;

%% Compute the global stuffs once
ComputeGlobal1D;

% flag to compute the misfit alone or with regularization
misfitFlag = 1; % misfit only
                % otherwise misfit + regularization

% Compute the misfit and its derivative at m
% [misfit, grad] = ThermalFinCost1Dmisfit(GamMap, misfitFlag);

% Sampling the gradient around the GamMap
m = length(GamMap);
N = 2000;
X = randn(m,N);
grads = zeros(m,N);
misfits = zeros(1,N);
gamma = 0.2;
for i=1:N
    [misfits(i), grads(:,i)] = ThermalFinCost1Dmisfit(GamMap+gamma*X(:,i), misfitFlag);
end

%% Find active subspace
[U,Sig,V] = svd(grads','econ');
sig = (1/m)*diag(Sig).^2;

%% Plots!
close all;
k=5;
figure;
plot(X'*V(:,1),U(:,k),'bx');
set(gca,'FontSize',14,'LineWidth',2);
axis square; grid on;

figure;
scatter(X'*V(:,1),X'*V(:,2),150,U(:,k),'filled');
axis square;



%%
close all;
k=30;
a = [ones(N,1) X']\U(:,k);
u = a(2:end)/norm(a(2:end));
figure;
plot(X'*u,U(:,k),'bx');
set(gca,'FontSize',14,'LineWidth',2);
axis square; grid on;









