function mu = ComponentViscosity(T, type, coeffs)

% ComponentViscosity - Calculates the viscosity of each species
%
% Syntax:  mu = ComponentViscosity(T, type, coeffs)
%
% Inputs:
%    T - temperature (in K)
%    type - type of correlation for calculating the vapor pressure
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    mu - viscosity (in cP)
%
% --------------------------- BEGIN CODE -------------------------------- %

    mu = zeros(1,length(type));
    for i=1:length(type)
        if (strcmp(type(i), 'law-1'))
            mu(i) = 10^(coeffs{i}(1)+coeffs{i}(2)/T+coeffs{i}(3)*T+coeffs{i}(4)*T^2);
        else
            error('Only the type=law-1 is currently available for viscosity');
        end
    end 

end
