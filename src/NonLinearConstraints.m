function [c, ceq] = NonLinearConstraints(z)

% NonLinearConstraints - Non linear constraints for surrogate optimization
% Not every MATLAB optimizer can handle non linear constraints. In the
% current version only 'ga' (genetic algorithm) manages non linear
% constraints
%
% Syntax:  [c, ceq] = NonLinearConstraints(z)
%
% Inputs:
%    z - Current optimization variables
%
% Outputs:
%    c - vector of non linear inequalities
%    ceq - vector of non linear equalities
%
% --------------------------- BEGIN CODE -------------------------------- %

    global target;
    global database;
    global optimizer;
    
    % Number of components
    nx = database.nx;

    % Integer variable: 0=species is not present, 1=species is present
    family = z(1:nx);
    
    % Molar fractions
    x = z((nx+1):2*nx).*family;
    sumx = sum(x);
    x = x/sumx;
    
    % Molecular Weight
    MW_surrogate = sum(x.*database.nH)*1.+sum(x.*database.nC)*12.;

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
    
    % Minimum mole fraction
    
    
    
    
    % Non-linear constraints
    c = [];
    
    % 1. Minimum amount of a single species
    if (target.minx ~= 0.)
        
        min_x = min(x(x ~= 0));
        
        c = [ -min_x+target.minx; ];
    
    end
    
    % 1. Volumetric fractions for families
    if (target.v_constraints == true)

        phi_families = database.family_table*V';
        
        c_v = zeros(length(phi_families)*2,1);
        for i=1:length(phi_families)
            c_v(2*i-1) =  phi_families(i)-target.uv(i);
            c_v(2*i)   = -phi_families(i)+target.lv(i);
        end
        
        c = [c; c_v];
        
    end

    
    % 3. Relative range interval for MW, HC, CN, and TIS 
    if (optimizer.rel_constraints == true)
        
        c_rel = [    MW_surrogate/target.MW-(1.+optimizer.delta_rel_MW);
                    -MW_surrogate/target.MW+(1.-optimizer.delta_rel_MW);
                     HC_surrogate/target.HC-(1.+optimizer.delta_rel_HC);
                    -HC_surrogate/target.HC+(1.-optimizer.delta_rel_HC);  
                     CN_surrogate/target.CN-(1.+optimizer.delta_rel_CN);
                    -CN_surrogate/target.CN+(1.-optimizer.delta_rel_CN);
                     TSI_surrogate/target.TSI-(1.+optimizer.delta_rel_TSI);
                    -TSI_surrogate/target.TSI+(1.-optimizer.delta_rel_TSI); 
                ];

        c = [c; c_rel];
        
    end
        
    % 4. Absolute range interval for MW, HC, CN, and TIS   
    if (optimizer.abs_constraints == true)
        
        c_abs = [    MW_surrogate-(target.MW+optimizer.delta_abs_MW);
                    -MW_surrogate+(target.MW-optimizer.delta_abs_MW); 
                     HC_surrogate-(target.HC+optimizer.delta_abs_HC);
                    -HC_surrogate+(target.HC-optimizer.delta_abs_HC);   
                     CN_surrogate-(target.CN+optimizer.delta_abs_CN);
                    -CN_surrogate+(target.CN-optimizer.delta_abs_CN); 
                     TSI_surrogate-(target.TSI+optimizer.delta_abs_TSI);
                    -TSI_surrogate+(target.TSI-optimizer.delta_abs_TSI);  
                ];
            
        c = [c; c_abs];
        
    end
              
    ceq = [];

end