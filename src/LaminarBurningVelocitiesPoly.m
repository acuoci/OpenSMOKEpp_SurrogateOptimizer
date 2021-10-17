function lbv = LaminarBurningVelocitiesPoly(lbvpures, Coeffs, exponents, x0)
% LaminarBurningVelocities - Calculates the LBV curve using a mixing rule 
%                            that is of polynomial type 
%
% Syntax:  lbv = ... 
%          LaminarBurningVelocitiesPoly(lbvpures, Coeffs, exponents, x0)
%
% Inputs:
%    v_pures - laminar flame speed for the pure components
%    Coeffs - optimized coefficients for the mixing rule
%    exponents - matrix of exponents (assign coefficients with their power)
%    x0 - mixture composition (in molar fractions)
%
% Output:
%    lbv - vector of ignition delay times (cm/s)
%
% --------------------------- BEGIN CODE -------------------------------- %

 lbv = zeros(size(lbvpures,1),1);
for a = 1:size(exponents,1)
    piece = ones(size(lbvpures,1),1);
    for b = 1:size(exponents,2)
        piece = piece.* (x0(b).*lbvpures(:,b)).^exponents(a,b);
    end
    lbv = lbv + Coeffs(a).*piece;
end
    

end
