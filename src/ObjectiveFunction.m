function f = ObjectiveFunction(z)

    global database;
    global target;
    global optimizer;
    
    % Number of components
    nx = database.nx;

    if (strcmp(optimizer.type,'ga'))
        
        % Integer variable: 0=species is not present, 1=species is present
        family = z(1:nx);

        % Molar fractions
        x = z((nx+1):2*nx).*family;
        sumx = sum(x);
        x = x/sumx;
        
    else
        
        x = z;
        sumx = sum(x);
        x = x/sumx;
        
    end   
    
    % Molecular Weight
    MW_surrogate = sum(x.*database.nH)*1+sum(x.*database.nC)*12;
    
    % H/C ratio
    HC_surrogate = sum(x.*database.nH)/sum(x.*database.nC);
    
    % Mass fractions
    omega = x.*database.MW/MW_surrogate;
    
    % Density (kg/m3)
    rhos = ComponentDensity(database.rhoType, database.rhoCoeffs);
    
    % Volumetric fractions (kg/m3)
    V = omega./rhos;
    sumV = sum(V);
    V = V/sumV;
    
    % Cetane number
    CN_surrogate = V*database.CN';
    
    % Cetane number
    TSI_surrogate = x*database.TSI';
    
    % Distillation curve
    optimizer.F_DC = 0.;
    if (optimizer.weight_DC ~= 0)
        
        % Calculates the distillation curve
        [vol, Td] = DistillationCurve(target.P, x, target.delta, rhos, database.MW, database.vpType, database.vpCoeffs);
        Tdi_target = interp1q(target.vol, target.Td, target.voli');
        Tdi_current = interp1q(vol', Td', target.voli');
        
        % Calculates the error (see Narayanaswamy and Pepiot, Comb. Theory Mod. 5, p. 883-897, 2018)
        sigma2 = target.sigma_DC^2;
        for i=1:length(target.voli)
            optimizer.F_DC = optimizer.F_DC + (1.-exp(-(Tdi_target(i)-Tdi_current(i))^2./2/sigma2))^2.;
        end
        optimizer.F_DC = optimizer.F_DC/length(target.voli);
        
    end
    
    % Objective function (single contributions)
    if (optimizer.error_type == 1)
        
        optimizer.F_HC  = (1. - exp( -(HC_surrogate-target.HC)^2/(2*target.sigma_HC^2)))^2;
        optimizer.F_MW  = (1. - exp( -(MW_surrogate-target.MW)^2/(2*target.sigma_MW^2)))^2;
        optimizer.F_CN  = (1. - exp( -(CN_surrogate-target.CN)^2/(2*target.sigma_CN^2)))^2;
        optimizer.F_TSI  = (1. - exp( -(TSI_surrogate-target.TSI)^2/(2*target.sigma_TSI^2)))^2;
        
    elseif (optimizer.error_type == 2)
        
        optimizer.F_HC  = (1. - HC_surrogate/target.HC)^2;
        optimizer.F_MW  = (1. - MW_surrogate/target.MW)^2;
        optimizer.F_CN  = (1. - CN_surrogate/target.CN)^2;
        optimizer.F_TSI = (1. - TSI_surrogate/target.TSI)^2; 
        
    elseif (optimizer.error_type == 3)
        
        optimizer.F_HC  = (1. - exp( -(HC_surrogate-target.HC)^2/target.HC))^2;
        optimizer.F_MW  = (1. - exp( -(MW_surrogate-target.MW)^2/target.MW))^2;
        optimizer.F_CN  = (1. - exp( -(CN_surrogate-target.CN)^2/target.CN))^2;
        optimizer.F_TSI = (1. - exp( -(TSI_surrogate-target.TSI)^2/target.TSI))^2;
         
    end

    % Objective function
    f = optimizer.weight_HC*optimizer.F_HC + ...
        optimizer.weight_MW*optimizer.F_MW + ...
        optimizer.weight_CN*optimizer.F_CN + ...
        optimizer.weight_TSI*optimizer.F_TSI + ...
        optimizer.weight_DC*optimizer.F_DC;

end