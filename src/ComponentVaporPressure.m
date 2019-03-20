% Calculate and return vapor pressure from parameters Antoine and temperature
function vp = ComponentVaporPressure(T, type, coeffs)

    vp = zeros(1,length(type));
    for i=1:length(type)
        if (strcmp(type(i), 'Antoine3'))
            vp(i) = 10^(coeffs{i}(1)-coeffs{i}(2)/(coeffs{i}(3)+T));
        else
            error('Only the type=Antoine3 is currently available for vapor pressure');
        end
    end 

end