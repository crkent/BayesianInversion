% Number of shape functions for solution interpolation
nshapeSolution = GetNoShapesLine(SolutionOrder);

% Number of elements on each region
nelem = size(mesh.ElementGroup,1);

% Number of Nodes for each element
nelemNode = size(mesh.ElementGroup,2);

% Finding the size of the stiffness matrix
nsize = mesh.nodes;

% Gauss point
[sg, ng, wg] = QuadLine(SolutionOrder);

Jacobian = zeros(ng,nelem);

DerivativeShape=zeros(nshapeSolution,ng,nelem);

Shape = zeros(nshapeSolution, ng);

Node = zeros(nelemNode,1);
x    = zeros(nelemNode,1);

for ie = 1:nelem
  % Take global node for current element
  Node(1:nelemNode) = mesh.ElementGroup(ie, 1:nelemNode);

  % x-y coordinates for each global node
  x(1:nelemNode) = mesh.coor(Node(1:nelemNode));
  
  for ig=1:ng
    [Shape(:,ig), DerivativeShape(:,ig,ie), Jacobian(ig,ie)]= shape_line(SolutionOrder, GeometryOrder, sg(:,ig), x);
    
  end
end

