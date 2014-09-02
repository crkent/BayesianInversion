function [sg, ng, wg]=QuadLine(order,nsflag)

% Finding the number of Gauss point corresponding to the order

% order    : The order of interpolation
% nsflag   : 2 means \int{\phi_i * \phi_j}
%            1 means \inf{\phi_i}

if 1,%(nsflag == 2)
    switch order
% $$$         case 1,
% $$$             ng = 2;
% $$$             sg = [-1.0/sqrt(3.), 1./sqrt(3.)];
% $$$             wg = [1., 1.];
% $$$      case 1,
% $$$       ng = 1;
% $$$       sg = [0.0];
% $$$       wg = [2.0];
     case {1,2},
      ng = 3;
      sg = [-sqrt(3./5.),  0.0,  sqrt(3./5.)];
      wg = [5./9.,  8./9., 5./9.];
     case {3},
      ng = 4;
      sg = [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7),...
           sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)];
      wg = [(18 - sqrt(30))/36, (18 + sqrt(30))/36,...
            (18 + sqrt(30))/36, (18 - sqrt(30))/36];
     otherwise
      disp('I am not supporting this case');
      disp('Please add order quadrature rules.')
      termination = stop
    end
else
    if (nsflag == 1)
        switch order
            case 1,
                ng = 1;
                sg = [0.0];
                wg = [2.0];
            case 2,
                ng = 2;
                sg = [-1.0/sqrt(3.), 1./sqrt(3.)];
                wg = [1.0, 1.0];
            otherwise
                disp('I am not supporting this case');
                disp('Please add order quadrature rules.')
                termination = stop
        end
    else
        disp('I am not supporting this case');
        disp('Please add boundary terms here.')
        termination = stop
    end
end
    

