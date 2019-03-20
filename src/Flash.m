function y=Flash(T, P, x0, alpha, type, coeffs)

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