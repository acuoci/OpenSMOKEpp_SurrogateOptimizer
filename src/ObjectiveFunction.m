function [f,F] = ObjectiveFunction(z,allowedspecies,database,target,optimizer,net,idtkriging,lbvmodel)

% ObjectiveFunction - Calculates the objective function. The optimization
% variables z are the number of moles of individual components. In case a
% a genetic algorithm is adopted, the vector of variables z is divided in
% two parts: the first part corresponds to integers (0 or 1) indicating if
% a species is present, while the second part corresponds to the number of
% moles of each species
%
% Syntax:  f = ObjectiveFunction(z,allowedspecies,...)
%
% Inputs:
%    z - Current optimization variables
%    allowedspecies - A logical vector that represents the species that can
%                     be used (1) or are excluded by default (0)
% Outputs:
%    f - objective function
%    F - a struct with the contributions and weights of each property to
%        the objective function
% --------------------------- BEGIN CODE -------------------------------- %

    % Number of components
    nx = database.nx;

    if (strcmp(optimizer.type,'ga'))
        
        % Integer variable: 0=species is not present, 1=species is present
        family = z(1:nx);

        % Molar fractions
        x = z((nx+1):2*nx).*family;
        x = x/sum(x);
        if sum(family)==0, x=zeros(1,nx); end
        xtemp = x;
    else
        x = zeros(1,length(allowedspecies));
        xtemp = z;
    end
    iz = 1;
    for ix = 1:length(allowedspecies)
        if allowedspecies(ix) == 1
            x(ix) = xtemp(iz);
            iz = iz+1;
        end
    end
    x = x/sum(x);  
    
    % Molecular Weight
    MW_surrogate = sum(x.*database.nH)*1+sum(x.*database.nC)*12;
    
    % H/C ratio
    HC_surrogate = sum(x.*database.nH)/sum(x.*database.nC);
    
    % Mass fractions
    omega = x.*database.MW/MW_surrogate;
    
    % Density (kg/m3)
    rhos = ComponentDensity(database.rhoType, database.rhoCoeffs);
    rho_surrogate=x*rhos';
    
    % Volumetric fractions (kg/m3)
    V = omega./rhos;
    sumV = sum(V);
    V = V/sumV;
    
    % Cetane number
    CN_surrogate = V*database.CN';
    
    % TSI number
    TSI_surrogate = x*database.TSI';
    
    % Viscosity
    mus = ComponentViscosity(target.Tvis, database.muType, database.muCoeffs);    % Viscosity single components
    mu_surrogate = sum(omega.*(mus.^(1/3)))^3;                                    % Blending rule: mu^1/3 = sum (omega_i*mu_i^1/3)
    
    % YSI
    YSI_surrogate=omega*database.YSI';
    
    % Distillation curve
    F_DC = 0.;
    if (optimizer.weight_DC ~= 0)
        
        % Calculates the distillation curve
        [vol, Td] = DistillationCurve(target.P, x, target.delta, rhos, database.MW, database.vpType, database.vpCoeffs);
        Tdi_target = interp1(target.vol, target.Td, target.voli');
        Tdi_current = interp1(vol, Td, target.voli');
        
        % Calculates the error 
        for i=1:length(target.voli)
            %F_DC = F_DC + (1.-exp(-(Tdi_target(i)-Tdi_current(i))^2./2/sigma2))^2.;
            F_DC = F_DC + (1-Tdi_current(i)/Tdi_target(i))^2;
        end
       F_DC = F_DC/length(target.voli);
        
    end

    % Ignition delay times
    if optimizer.weight_idt==0
        idt_surrogate = 0;
    else
        if target.IdtCalculationType == "kriging"
        idt_surrogate = IgnitionDelayTimes(target.list_species,target.TIdt, target.PIdt, target.phiIdt,...
         x, idtkriging.species, idtkriging.model, target.tauType, target.IdtCalculationType, idtkriging.logarithm); 
        else
        idt_surrogate = IgnitionDelayTimes(target.list_species,target.TIdt, target.PIdt, target.phiIdt,...
         x, net.species, net.net, target.tauType, target.IdtCalculationType, net.logarithm);
        end
    end
    F_IDT = 0;
    
    % Laminar burning velocities
    if optimizer.weight_lbv==0
        lbv_surrogate = 0;
    else
        LBVpures = lbvmodel.sLpures;
        LBVtarget = target.lbv;
        try
        if numel(lbvmodel.phipures) ~= numel(target.phiLBV)
            % if phi do not coincide, do an interpolation
            phimin = max([min(lbvmodel.phipures),min(target.phiLBV)]);
            phimax = min([max(lbvmodel.phipures),max(target.phiLBV)]);
            idxphi = find( (lbvmodel.phipures >= phimin) .* (lbvmodel.phipures <= phimax) );
            phis = lbvmodel.phipures(idxphi);
            LBVpures = lbvmodel.sLpures(idxphi,:);
            LBVtarget = interp1(target.phiLBV,target.lbv,phis);
        end
        catch
        end
        lbv_surrogate = LaminarBurningVelocitiesPoly(LBVpures, ...
               lbvmodel.LBVcoeffs, lbvmodel.LBVexponents, x);
%         lbv_surrogate = LaminarBurningVelocitiesPoly(lbvmodel.sLpures, ...
%                lbvmodel.LBVcoeffs, lbvmodel.LBVexponents, x);
    end
    F_LBV = 0;
    
    % Objective function (single contributions)
        F_HC  = (1. - HC_surrogate/target.HC)^2;
        F_MW  = (1. - MW_surrogate/target.MW)^2;
        F_CN  = (1. - CN_surrogate/target.CN)^2;
        F_TSI = (1. - TSI_surrogate/target.TSI)^2; 
        F_mu  = (1. - mu_surrogate/target.mu)^2;
        F_YSI = (1. - YSI_surrogate/target.YSI)^2;
        F_rho = (1. - rho_surrogate/target.rho)^2;
        if optimizer.weight_idt~=0
        for i=1:length(idt_surrogate)
            F_IDT = F_IDT + (1. - log(idt_surrogate(i))/log(target.Idt(i)))^2; 
        end
        F_IDT = F_IDT./numel(target.Idt);
        end
        if optimizer.weight_lbv~=0
        for i=1:length(lbv_surrogate)
            F_LBV = F_LBV + (1. - lbv_surrogate(i)/LBVtarget(i))^2;
        end
        F_LBV = F_LBV./numel(target.lbv);
        end

    F.F_HC  = F_HC;    F.F_MW  = F_MW;    
    F.F_CN  = F_CN;    F.F_TSI = F_TSI;
    F.F_mu  = F_mu;    F.F_YSI = F_YSI;   
    F.F_rho = F_rho;   F.F_IDT = F_IDT;
    F.F_LBV = F_LBV;   F.F_DC  = F_DC;
    F.w_HC  = optimizer.weight_HC;    F.w_MW  = optimizer.weight_MW;
    F.w_CN  = optimizer.weight_CN;    F.w_TSI = optimizer.weight_TSI;
    F.w_mu  = optimizer.weight_mu;    F.w_YSI = optimizer.weight_YSI;
    F.w_rho = optimizer.weight_rho;   F.w_IDT = optimizer.weight_idt;
    F.w_LBV = optimizer.weight_lbv;   F.w_DC  = optimizer.weight_DC;
    % Objective function
    f = optimizer.weight_HC*F_HC + ...
        optimizer.weight_MW*F_MW + ...
        optimizer.weight_CN*F_CN + ...
        optimizer.weight_TSI*F_TSI + ...
        optimizer.weight_mu*F_mu + ...
        optimizer.weight_YSI*F_YSI + ...        
        optimizer.weight_DC*F_DC + ...
        optimizer.weight_rho*F_rho + ...
        optimizer.weight_idt*F_IDT + ...
        optimizer.weight_lbv*F_LBV;

end