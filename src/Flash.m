function T=Flash(Tguess, P, x0, alpha, type, coeffs)
% Flash - Solves a flash
%
% Syntax:  y=Flash(T, P, x0, alpha, type, coeffs)
%
% Inputs:
%    T - first guess temperature (in C)
%    P - pressure (in mmHg)
%    x0 - mixture composition (in molar fractions)
%    alpha - 
%    type - type of correlation for calculating the vapor pressure
%    coeffs - vector of coefficients from the database
%
% Outputs:
%    T - temperature
%
% --------------------------- BEGIN CODE -------------------------------- %

%     options = optimset('Display','off');
%     T = fsolve(@balance, Tguess, options);
    options.tolX = 1e-4; 
    options.tolF = 1e-6; 
    options.maxiter = 250;
    [T,it] =  AzzeraBalanceSecanti(Tguess-100,Tguess+100,type,coeffs,P,x0,alpha,options);
    if it == options.maxiter
        options.maxiter = 500;
        [T,it] =  AzzeraBalanceSecanti(90,230,type,coeffs,P,x0,alpha,options);
        if it == options.maxiter
			fprintf('Debug: flash done again, and it''s not enough\n');
        else
%             fprintf('O');
%             fprintf('Debug: flash done again (250+%d)\n',it);
        end
        % If it does not converge at first time, enlarge the range of
        % search between 90°C and 130°C - the max and min boiling points of
        % pure components - and allow more iterations
    else
%         fprintf('Debug: flash done once (%d)\n',it);
%         fprintf('-');
    end

    % Equation to be zeroed to close the mole balance 
    % temperature, liquid molar composition, fraction that vaporized and pressure
    function f = balance(Temp,type,coeffs,P,x0,alpha)             

            K = ComponentVaporPressure(Temp, type, coeffs)/P;
            f=x0.*(K-1.)./(1.+alpha*(K-1.));
            f=sum(f);

    end
	
    function [T2,iter] = AzzeraBalanceSecanti(T1,T2,type,coeffs,P,x0,alpha,options)
        y1 = balance(T1,type,coeffs,P,x0,alpha);
        y2 = balance(T2,type,coeffs,P,x0,alpha);

        m = (y2-y1)/(T2-T1);
        T1 = T2;
        T2 = T2 - y2/m;
        y1 = y2;
        for iter = 1:options.maxiter
            y2 = balance(T2,type,coeffs,P,x0,alpha);

            if abs(y2)<options.tolF || abs(T1-T2)<options.tolX
                T = T2;
                break;
            end
            m = (y2-y1)/(T2-T1);
            T1 = T2;
            T2 = T2 - y2/m;
            y1 = y2;
        end
    end

end
