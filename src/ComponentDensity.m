function rho = ComponentDensity(type, coeffs)

    rho = zeros(1,length(type));
    for i=1:length(type)
        if (strcmp(type(i), 'constant'))
            rho(i) = coeffs(1,i);
        else
            error('Only the type=constant is currently available for density');
        end
    end