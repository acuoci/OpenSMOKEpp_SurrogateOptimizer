% Calculate and return 
function Tb = ComponentBoilingTemperature(p, type, coeffs)

% ComponentBoilingTemperature - Calculates the boiling temperature
% from Antoine parameters and temperature
%
% Syntax:  Tb = ComponentBoilingTemperature(p, type, coeffs)
%
% Inputs:
%    p - pressure (in mmHg)
%    type - type of correlation for calculating the vapor pressure
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    Tb - vapor pressure (in mmHg)
%
% --------------------------- BEGIN CODE -------------------------------- %

    Tb = zeros(1,length(type));
    for i=1:length(type)
        if (strcmp(type(i), 'Antoine3'))
            Tb(i) = coeffs{i}(2)/(coeffs{i}(1)-log10(p)) - coeffs{i}(3);
        else
            error('Only the type=Antoine3 is currently available for vapor pressure');
        end
    end 

end