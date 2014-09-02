function [shape, dnds] = ShapeDerivativeLine(nshape,s);

% nshape     : number of shape functions
% s          : Gauss point

% shape      : Values of shape functions
% dnds       : derivative of shape function wrt s_i

switch nshape
    case 1,
        disp('I am not supporting this case');
        disp('Constant geometry and solution interpolation are not used.')
        termination = stop
    case 2,
        %     shp(3,j): value of shape function j at quadrature points
        shape(1,1) = 0.5*(1.-s);
        shape(2,1) = 0.5*(1.0 + s);

        %     The derivatives of shape function wrt reference coors dN/ds
        %     dnds(i,j): derivative of shape function j wrt s_i

        dnds(1,1) = -0.5;
        
        dnds(2,1) = 0.5;
        
    case 3,
        %     shp(3,j): value of shape function j at quadrature points
        shape(1,1) = 0.5*(1.-s) - 0.5*(1-s)*(1+s);
        shape(2,1) = 0.5*(1.+ s) - 0.5*(1-s)*(1+s);
        shape(3,1) = (1-s)*(1+s);

        dnds(1,1) = s - 0.5;

        dnds(2,1) = s + 0.5;
        
        dnds(3,1) = -2.*s;
        
    otherwise,
        disp('I am not supporting this case');
        disp('You can add higher interpolation here.')
        termination = stop
end