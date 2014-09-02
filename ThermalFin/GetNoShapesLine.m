function [nshape] = GetNoShapesLine(order)

% order    : Solution order on each elemment

switch order
    case 0,
        nshape = 1;
    case 1,
        nshape = 2;
    case 2,
        nshape = 3;
    otherwise,
        disp('I am not supporting this case');
        disp('Please add order solution order.')
        termination = stop
end