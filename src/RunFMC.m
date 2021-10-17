function solution = runFMC(allowedspecies,x0,optimizer,target,f_min_con,database,net,idtkriging,lbvmodel,isimu,totsimu)
   if nargin < 10, isimu = 0; totsimu = 0; end
    optimizer.type='fmincon';
    optimizer.opts = optimoptions(  'fmincon', ...
        'MaxIterations', 5000, ...
        'Algorithm', 'interior-point', ... 'interior-point' | 'trust-region-reflective' | 'sqp' | 'active-set'
        'Display','off',...
        'StepTolerance', 1e-12,...
        'ConstraintTolerance',1e-12);
        
    %optimizer.opts.UseParallel=optimizer.parallel_computing;
    optimizer.n = database.nx;
    optimizer.lb = zeros(1,nnz(allowedspecies));
    optimizer.ub = optimizer.lb + 1;
    a = find(allowedspecies);
    xFG = x0(a);
    
%     % Add random noise
%     xFG = xFG + randn(size(xFG)).*0.1;
%     xFG = xFG/sum(xFG);
    
    % Optimize
    try
    [sol,fval,exitflag,output] = fmincon(@(z)ObjectiveFunction(z,allowedspecies,database,target,optimizer,net,idtkriging,lbvmodel), ...
                                        xFG, [], [], [], [], optimizer.lb, optimizer.ub, [], ...
                                        optimizer.opts);
    solution.z = sol;
    solution.x = zeros(1,length(allowedspecies));
    iz = 1;
    for ix = 1:length(allowedspecies)
        if allowedspecies(ix) == 1
            solution.x(ix) = sol(iz);
            iz = iz+1;
        end
    end
    solution.x = solution.x./sum(solution.x);
    % Remove x below 1/1000
    toremove = (solution.x < 0.001);
    solution.x(toremove) = 0;
    solution.x = solution.x./sum(solution.x);
    solution.fval = fval;
    solution.exitflag = exitflag;
    solution.output = output;
    solution.done = true;
    catch ex
    solution = [];
    solution.done = false;
    solution.output = ex;
    solution.exitflag = -100;
    solution.x = zeros(1,length(allowedspecies)) + NaN;
%     fprintf('Some error occured: %s\n',ex.message);
    end
    solution.solvername = 'fmincon';
    solution.allowedspecies = allowedspecies;

    fprintf('  [%3.3d/%3.3d]\tfmincon       end - flag %d\n',isimu,totsimu,solution.exitflag);
end