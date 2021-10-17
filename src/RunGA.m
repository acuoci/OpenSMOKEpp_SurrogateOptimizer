function solution = runGA(optimizer,target,genetic_algorithm,database)
    global net idtkriging lbvmodel
    % Optimizer variables
    optimizer.n = 2*database.nx;            % total number of variables
    optimizer.ub = ones(1,optimizer.n);     % upper boundaries
    optimizer.ub(database.nx+1:end) = Inf;  % upper boundaries
    optimizer.lb = zeros(1,optimizer.n);    % lower boundaries
    optimizer.intcon = 1:1:database.nx;     % integer unknowns
    
    % Options
    optimizer.opts = optimoptions('ga',     'MaxStallGenerations',genetic_algorithm.maxStallGen,...
                                            'FunctionTolerance',genetic_algorithm.FunTol,...
                                            'PopulationSize', genetic_algorithm.population_size, ...
                                            'MaxGenerations',genetic_algorithm.maxGen,...
                                            'Display','off',...
                                            ...'PlotFcn',{@gaplotbestf,@gaplotstopping},...
                                            'CrossoverFraction',genetic_algorithm.crossover_fraction);
    if optimizer.parallel_computing == true
        optimizer.opts.UseParallel=true;
        optimizer.opts.UseVectorized;
    end

    % Check consistency of input variables
    target.nmin = sum(target.ln);
    if (target.nmax < target.nmin)
        error('Error in target definition: inconsistency between min and max number of components');
    end

    % Check consistency of input variables
    target.vmax = sum(target.uv);
    if (target.vmax <= 1)
        error('Error in target definition: inconsistency between upper boundaries for volumetric fractions');
    end

    % Linear constrants on number of species
    A = [];
    b = [];
    if (target.global_constraints_nspecies == true)
        A = [  A;
               ones(1,database.nx),  zeros(1,database.nx); ...
              -ones(1,database.nx), -zeros(1,database.nx);  ];
        b = [  b;
               target.nmax; ...
              -target.nmin; ];
    end
    if (target.family_constraints_nspecies == true)
        A = [  A;
               [database.family_table'; zeros(database.nx,length(family_names))]'; 
              -[database.family_table'; zeros(database.nx,length(family_names))]'; ];
        b = [  b;
               target.un';
              -target.ln';];
    end

    % Optimize
    [sol,fval,exitflag,output] =    ga(@(z)ObjectiveFunction(z,target.active_species,database,target,optimizer,net,idtkriging,lbvmodel),...
                                    optimizer.n, A, b, [],[], ...
                                    optimizer.lb,optimizer.ub, ...
                                    @NonLinearConstraints,optimizer.intcon, ...
                                    optimizer.opts);

    % Solution
    solution.n = sol(1:database.nx).*sol((database.nx+1):2*database.nx);
    sumn = sum(solution.n);
    solution.x = solution.n/sumn;
    solution.fval = fval;
    solution.exitflag = exitflag;
    solution.output = output;
    solution.solvername = 'ga';
end