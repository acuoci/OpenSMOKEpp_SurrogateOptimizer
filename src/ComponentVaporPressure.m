% Calculate and return 
function vp = ComponentVaporPressure(T, type, coeffs)

% ComponentVaporPressure - Calculates the vapor pressure from Antoine 
% parameters and temperature
%
% Syntax:  vp = ComponentVaporPressure(T, type, coeffs)
%
% Inputs:
%    T - temperature (in C)
%    type - type of correlation for calculating the vapor pressure
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    rho - vapor pressure (in mmHg)
%
% --------------------------- BEGIN CODE -------------------------------- %

    vp = zeros(1,length(type));
    for i=1:length(type)
        if (strcmp(type(i), 'Antoine3'))
            vp(i) = 10^(coeffs{i}(1)-coeffs{i}(2)/(coeffs{i}(3)+T));
        else
            error('Only the type=Antoine3 is currently available for vapor pressure');
        end
    end 

end