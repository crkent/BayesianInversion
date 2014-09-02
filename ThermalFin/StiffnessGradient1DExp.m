function [G] = StiffnessGradient1DExp(Gam);
% This function will take the configuration and a triangulation and return Aq and Fh
% by using finite element method

global nelem nelemNode
global nshapeSolution
global mesh data
global wg Jacobian DerivativeShape Shape

G = zeros(size(Gam));
    
Node = zeros(nelemNode,1);
   
% Loop over all elements
for ie = 1:nelem

  % Take global node for current element
  Node(1:nelemNode) = mesh.ElementGroup(ie, 1:nelemNode);
  
  %% Compute the element matrix using quadrature rule
  %% This is actually the Galerkin term
  gam = Gam(Node);
  U = data.ForwardSolution(Node);
  V = data.AdjointSolution(Node);
%  go =  ElementStiffnessGradient1DExp(ie, gam, U, V);

  %----------- Compute the gradient-----------------------
  conductivity =  exp(Shape.' * gam);
  
  % compute the Gauss weights
  w = wg.'.*Jacobian(:,ie);
  
  Derivative = DerivativeShape(:,:,ie);
  
  GradU_GradV = w.* conductivity .* ((Derivative.'*U) .* (Derivative.'*V));
  g = Shape * GradU_GradV;
  %------------------------------------------------------ 
  
  % ensemble the gradient
  G(Node) = G(Node) + g;
% $$$   for alpha = 1:nshapeSolution
% $$$     i = Node(alpha);
% $$$     G(i) = G(i) + g(alpha);
% $$$   end
  
end