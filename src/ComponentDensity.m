function rho = ComponentDensity(type, coeffs)

% ComponentDensity - Calculates the density of each species
%
% Syntax:  rho = ComponentDensity(type, coeffs)
%
% Inputs:
%    type - type of correlation for calculating the density
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    rho - density (in kg/m3)
%
% --------------------------- BEGIN CODE -------------------------------- %

    rho = zeros(1,length(type));
    for i=1:length(type)
        if (strcmp(type(i), 'constant'))
            rho(i) = coeffs(1,i);
        else
            error('Only the type=constant is currently available for density');
        end
    end