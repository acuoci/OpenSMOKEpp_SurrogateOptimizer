function [vol, T_distillation] = ...
         DistillationCurve(P, x0, delta, rho, MW, type, coeffs)

% DistillationCurve - Calculates the distillation curve
%
% Syntax:  [vol, T_distillation] = ...
%          DistillationCurve(P, x0, delta, rho, MW, type, coeffs)
%
% Inputs:
%    P - pressure (in mmHg)
%    x0 - mixture composition (in molar fractions)
%    delta - step increment for evaluating the curve
%    rho - densities of species (in kg/m3)
%    MW - molecular weights of species (in kg/kmol)
%    type - type of correlation for calculating the vapor pressure
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    vol - vector of recovered volumes (%)
%    T_distillation - vector of distillation temperatures (in C)
%
% --------------------------- BEGIN CODE -------------------------------- %

    % Memory pre-allocation
    mol=0:delta:100;
    vol=0:delta:100;
    T_distillation=0:delta:100;

    % Initial mixture volume
    V0 = MoleFractions2MolarVolume(x0, rho, MW);
    
    % Bubble temperature
    T=0.;                                                   % Initial temperature (in C)
    T_distillation(1)=Flash(T, P, x0, 0., type, coeffs);    % First point in distillation curve

    % Loop over the whole range
    for i=1:length(mol)-1

        alpha=delta/(100.-mol(i));                          % Calculates the alpha to the next point
        T=Flash(T, P, x0, alpha, type, coeffs);             % Solves flash
        T_distillation(i+1)=T;                              % Stores T in distillation curve
        K = ComponentVaporPressure(T, type, coeffs)/P;      % Calculates K
        x0=x0./(1.+alpha*(K-1.));                           % Updates liquid composition
        y=K.*x0;                                            % Calculates the vapor composition 
        V = delta*MoleFractions2MolarVolume(y, rho, MW);    % Calculates the evaporated volume
        vol(i+1)=vol(i)+V/V0;                               % Updates the total evaporated volume 

    end

end

% Calculates the liquid volume
function Vtilde = MoleFractions2MolarVolume(x, rho, MW)

    Vtilde = sum(x.*MW./rho);
    
end

