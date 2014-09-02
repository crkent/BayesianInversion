function [Ah,Mh, Fh,Bh,BBh,Shape_obs] = FEM1DExp(Gam,mesh);
% This function will take the configuration and a triangulation and return Aq and Fh
% by using finite element method
global GeometryOrder SolutionOrder
global nshapeSolution nshapeSolutionLine
global nsize
global nelem nelemNode nseg nsegNode nsegRhs
global data
global wg Jacobian DerivativeShape Shape

%%%%%   INTERIOR ELEMENT CONTRIBUTION

ir = zeros(nsize * 2,1);
ic = zeros(nsize * 2,1);
s = zeros(nsize * 2,1);

% $$$ % Mass Matrix
% $$$ irm = zeros(nsize * 2,1);
% $$$ icm = zeros(nsize * 2,1);
% $$$ sm = zeros(nsize * 2,1);

Node = zeros(nelemNode,1);

ntriplet = 0;
% Loop over all elements
for ie = 1:nelem
    
  % Take global node for current element
  Node(1:nelemNode) = mesh.ElementGroup(ie, 1:nelemNode);
  
  %% Compute the element matrix using quadrature rule
  %% This is actually the Galerkin term
  gam = Gam(Node);
%  [AI,M]=ElementStiffness1D(ie,gam);
%  [AIo]=ElementStiffnessExp1D(ie,gam);
  
  %----------- Compute element stiffness--------------
  conductivity = exp(Shape.' * gam);

  w = wg.' .* Jacobian(:,ie);
  Derivative = DerivativeShape(:,:,ie);
  temp = diag(w .* conductivity);
  AI = (Derivative * temp) * Derivative.';

% $$$   if nargout > 1
% $$$     M  = (Shape * diag(w)) * Shape.';
% $$$   end
  %--------------------------------------------------
  
  % Assembly matrix A
  
  for alpha = 1:nshapeSolution
    i = Node(alpha);
    for beta = 1:nshapeSolution
      j = Node(beta);
      ntriplet = ntriplet + 1;
      len = length(ir);
      if ntriplet > len
        ir(2*len) = 0;
        rc(2*len) = 0;
        s(2*len) = 0;
        
% $$$         % Mass Matrix
% $$$         irm(2*len) = 0;
% $$$         rcm(2*len) = 0;
% $$$         sm(2*len) = 0;
      end
      ir(ntriplet) = i;
      ic(ntriplet) = j;
      s(ntriplet) = AI(alpha, beta);
      
% $$$       % Mass matrix
% $$$       irm(ntriplet) = i;
% $$$       icm(ntriplet) = j;
% $$$       sm(ntriplet) = M(alpha, beta);
    end
  end
end

if nargout > 1,
% $$$   % For testing the analytical solution
% $$$   Mh = sparse(irm(1:ntriplet),icm(1:ntriplet),sm(1:ntriplet),nsize, ...
% $$$             nsize);
  Mh = 0;
end

% contribution from the boundary
ntriplet = ntriplet + 1;
ir(ntriplet) = nsize;
ic(ntriplet) = nsize;
s(ntriplet)  = data.Bi;

Ah = sparse(ir(1:ntriplet),ic(1:ntriplet),s(1:ntriplet),nsize, ...
            nsize);

%%%%%%%%%     RHS COMPUTATION
if (nargout > 2)
  % Right hand side vector F
  Fh = zeros(nsize,1);    
  Fh(1) = 1;
end

if (nargout > 3)
  Shape_obs = zeros(nshapeSolution,data.Nobs);
  %%%%%%%%  Observation vector
  ntriplet = 0;
  irn = zeros(data.Nobs * nshapeSolution,1);
  icn = zeros(data.Nobs * nshapeSolution,1);
  sn = zeros(data.Nobs * nshapeSolution,1);
  % observation vector
  for ie = 1:data.Nobs
    % Take global node for current element
    Node(1:nelemNode) = mesh.ElementGroup(data.element_obs(ie), 1:nelemNode);
    [shapeSolution] = ShapeDerivativeLine(nshapeSolution,data.s_obs(ie));
    
    icn(nshapeSolution * (ie-1) + 1 : nshapeSolution * ie) = Node;
    irn(nshapeSolution * (ie-1) + 1 : nshapeSolution * ie) = ie;
    sn(nshapeSolution * (ie-1) + 1 : nshapeSolution * ie) = shapeSolution; 
    Shape_obs(:,ie) =  shapeSolution;
    for alpha = 1:nshapeSolution
      i = Node(alpha);
      for beta = 1:nshapeSolution
        j = Node(beta);
        ntriplet = ntriplet + 1;
        ir(ntriplet) = i;
        ic(ntriplet) = j;
        s(ntriplet) = shapeSolution(alpha) * shapeSolution(beta);
      end
    end
  end
  
  BBh = sparse(ir(1:ntriplet),ic(1:ntriplet),s(1:ntriplet),nsize, ...
               nsize);
  
  Bh = sparse(irn,icn,sn,data.Nobs, nsize);
  
end
