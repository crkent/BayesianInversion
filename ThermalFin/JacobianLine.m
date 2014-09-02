function [DetJacobian, InvJacobian] = JacobianLine(dndsGeometry,nshapeGeometry,x)

% dndsGeometry    : Derivatives of geometry interpolation functions (shapes)
% nshapeGeometry  : Number of shape function for geometry interpolation
% x               : x-y coordinates of the element

% DetJacobian     : The determinant of the Jacobian
% InvJacobian     : The inverse of the Jacobian


% Number of geometry nodes = Number of shape functions for geometry
nnode = size(x,1);

% dx/ds
dxds = x'*dndsGeometry;

% The Determinant of the Jacobian
DetJacobian = norm(dxds,2);

%% Fortran and C implementation

% %     compute the jacobian
% %     dxds(i,j): derivative of x_i wrt s_j
% for j = 1:1 % since we have 1 s_i
%     for  i = 1:2 % since we have 2 x_i
%         dxds(i,j) = 0.;
%         % transformation metric (2 x 3 matrix)
%         % since we have 3 s_i but x and y
%         for k = 1:nshapeGeometry
%             dxds(i,j) = dxds(i,j) + x(k,i)*dndsGeometry(k,j);
%         end
%     end
% end
% 
% % Notice that dS = sqrt(dx^2 + dy^2). But dx = (dx/ds)*ds and dy =
% % (dy/ds)*ds. Hence dS = sqrt((dx/ds)^2 + (dy/ds))*ds. As a result
% % The Jacobian is sqrt((dx/ds)^2 + (dy/ds)).
% 
% DetJacobian = sqrt(dxds(1,1)^2 + dxds(2,1)^2);

% Note that dx = (dx/ds)*ds, hence ds = (dx/ds)^{-1}*dx
% Hence (dx/ds)^{-1} is the inverse of Jacobian

%%% Let's not compute the Inverse of the transformation for now
%% Since the horizontal line or vertical line will have the singular
%%% tranformation

% InvJacobian(1,1) = 1.0/dxds(1,1);
% InvJacobian(1,2) = 1.0/dxds(2,1);

InvJacobian = 1/DetJacobian; %(not correct)