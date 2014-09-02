%%
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
[~,Sig,V] = svd(grads','econ');
evals = (1/m)*diag(Sig).^2;

%% Plots!
close all;
figure;
plot(X'*V(:,1),misfits,'bo','LineWidth',2);
set(gca,'FontSize',14);
axis square; grid on; xlabel('Active variable'); ylabel('Misfit');
%print('figs/misfits','-depsc','-r300');

figure;
plot(1:m,V(:,1),'rx',1:m,GamMap/norm(GamMap),'bo','LineWidth',2);
set(gca,'FontSize',14);
axis square; grid on; xlabel('Index'); ylabel('Direction');
%print('figs/evec','-depsc','-r300');

figure;
scatter(X'*V(:,1),X'*V(:,2),60,misfits,'filled');
axis square;
xlabel('y_1'); ylabel('y_2');

figure;
semilogy(1:m,evals,'rx','LineWidth',2);
set(gca,'FontSize',14); xlim([1 m]);
axis square; grid on; xlabel('Index'); ylabel('Eigenvalues');
%print('figs/eval','-depsc','-r300');









