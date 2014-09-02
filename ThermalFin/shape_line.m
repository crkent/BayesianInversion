function [Shape, Dshape, DetJacobian]=shape_line(SolutionOrder, GeometryOrder, s,x)  

% shape function for  3 node element
% SolutionOrder    : solution order of interpolation
% GeometryOrder    : geometry order of interpolation
% x:       : x-y coordinates of nodal points. Column k of x contains
%          : x-y coordinates of point k
% s        : Gausse points
% Jacobian : The Jacobian of the geometry transformation

% Number of geometry nodes
nnode = size(x,1);

% Number of shape functions for geometry interpolation
nshapeGeometry = GetNoShapesLine(GeometryOrder); 

% Number of shape functions for geometry interpolation <= nnode
if (nshapeGeometry > nnode)
    disp('Something wrong here.')
    disp('nshapeGeometry must be at most equal to nnode');
    disp('Please stop your program.')
    termination = stop
end

%Finding Shape derivatives for geometry
[shapeGeometry, dndsGeometry] = ShapeDerivativeLine(nshapeGeometry,s);

% Finding the Jacobian of the transformation
[DetJacobian, InvJacobian] = JacobianLine(dndsGeometry,nshapeGeometry,x(1:nshapeGeometry,:));

% Number of shape functions for solution interpolation
nshapeSolution = GetNoShapesLine(SolutionOrder); 

% Finding Shape derivatives for solution wrt s_i
[shapeSolution, dndsSolution] = ShapeDerivativeLine(nshapeSolution,s);

Shape = shapeSolution;

%%% Let's not compute the Inverse of the transformation for now
%% Since the horizontal line or vertical line will have the singular
%%% tranformation

% % Finding Shape derivatives for solution wrt x_i
% 
% %     The derivatives of shape function wrt the real coors dN/dx
% %     shp(i,j), i<3: derivative of shape function j wrt x_i
% for j = 1:nshapeSolution
%     for i = 1:2
%         Dshape(i,j) = 0.;
%         for k = 1:1 % we have 1 s_i
%             Dshape(i,j) = Dshape(i,j) + dndsSolution(k,j)*InvJacobian;%InvJacobian(k,i);
%         end
%     end
% end

if (size(x,2) == 2)
  Dshape = zeros(nshapeSolution,2);
elseif (size(x,2) == 1)
  Dshape = dndsSolution * InvJacobian;
else
  x
  error('not supported')
end