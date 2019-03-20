function y=Flash(T, P, x0, alpha, type, coeffs)

% Flash - Solves a flash
%
% Syntax:  y=Flash(T, P, x0, alpha, type, coeffs)
%
% Inputs:
%    T - temperature (in C)
%    P - pressure (in mmHg)
%    x0 - mixture composition (in molar fractions)
%    alpha - 
%    type - type of correlation for calculating the vapor pressure
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    y - molar composition
%
% --------------------------- BEGIN CODE -------------------------------- %

    options = optimset('Display','off');
    y = fsolve(@balance, T, options);

    % Equation to be zeroed to close the mole balance 
    % temperature, liquid molar composition, fraction that vaporized and pressure
    function f = balance(T)             

            K = ComponentVaporPressure(T, type, coeffs)/P;
            f=x0.*(K-1.)./(1.+alpha*(K-1.));
            f=sum(f);

    end

end