function [J,G] = ThermalFinCost1Dmisfit(Gam,misfitFlag)

global mesh data

if nargin < 2
  misfitFlag = 1;
end

%------Compute the misfit--------------------
% forward solve
[Ah] = FEM1DExp(Gam,mesh);
[L,U,P,Q] = lu(Ah);
% The exact QoI
uh = Q*(U\(L\(P * data.ForwardRHS)));
phi = uh(data.B);

% compute the cost function
J = 0.5 / data.sigma2 * sum((phi-data.phi_obs).^2);

% add regularization if needed
if misfitFlag ~= 1,
  prior = (data.PhiM * (Gam-data.Gam0)) ./ data.lam;
  J = J + 0.5*data.alpha * sum(prior.^2);
end
%--------------------------------------------

%--------compute the gradient of the misfit-----------
data.ForwardSolution = uh;
%  data.Stiffness = Ah;
% $$$ data.L = L; data.U = U; data.P = P; data.Q = Q;
data.phi = phi;

% Solve the adjoint problem
AdjointRHS = zeros(size(uh));
AdjointRHS(data.B) = -1/data.sigma2*(data.phi - data.phi_obs);
data.AdjointSolution =  Q*(U\(L\(P * AdjointRHS)));

G = StiffnessGradient1DExp(Gam);

if misfitFlag ~= 1,
  G = G + data.alpha  * data.PhiM.' * (prior ./ data.lam);
end  
%----------------------------------------------------

