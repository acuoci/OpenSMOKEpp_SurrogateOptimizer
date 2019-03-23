% ----------------------------------------------------------------------- %
%                   OPTIMIZER FOR FUEL SURROGATES                         %
% ----------------------------------------------------------------------- %
%                                                                         %
% This Matlab package implements an optimization algorithm to find the    %
% best fuel surrogate matching a set of experimental properties, such as  %
% molecular weight, H/C ratio, cetane number, TSI (Threshold Soot Index), %
% and distillation curves.                                                %
%                                                                         %
% The methodology is inspired by the following papers:                    %
%  * K. Narayanaswamy, H. Pitsch, P. Pepiot, Combustion and Flame,        %
%    165, p. 288-309 (2016)                                               %
%  * K. Narayanaswamy, P. Pepiot, Combustion Theory and Modelling,        %
%    22(5), p. 883-897 (2018)                                             %
%                                                                         %
% It is written purely in Matlab language. It is self-contained. There    %
% is no external dependency.                                              %
%                                                                         %
% Note: this package requires Matlab R2016b or latter, since it utilizes  %
% a new series of Matlab features.                                        %
%                                                                         %
% A. Cuoci, D. Demetryouss, M. Mehl                                       %
% CRECK Modeling Lab                                                      %
% Department of Chemistry, Materials, and Chemical Engineering            %
% Politecnico di Milano (Italy)                                           %
% ----------------------------------------------------------------------- %

clear variables;
close all;

global target;          % definition of targets for optimization
global database;        % database of individual components
global optimizer;       % optimizer options

% ----------------------------------------------------------------------- %
% Species available in the database
% ----------------------------------------------------------------------- %
% n-alkanes:  'NC4H10' 'NC5H12' 'NC7H16' 'NC10H22' 'NC12H26' 'NC14H30'
% i-alkanes:  'IC5H12' 'IC8H18'
% c1-alkanes: 'CHX', 'MCH'
% aromatics:  'T124MBZ', 'O-XYL'


% ----------------------------------------------------------------------- %
% [UD] Optimization Setup 
% ----------------------------------------------------------------------- %

% Optimization algorithm
% 1. 'ga': genetic algorithm
% 2. 'particleswarm': particle swarm solver (derivative-free)
% 3. 'patternsearch': pattern search solver (derivative-free)
optimizer.type = 'ga';

% Names of families
family_names = { 'n-alkanes', 'i-alkanes', 'c1-alkanes', 'aromatics' };

% List of active species
list_species = { 'NC4H10', 'NC10H22', 'NC16H34', ...
                 'IC5H12', 'IC8H18', ...
                 'CHX', 'MCH', ...
                 'T124MBZ', 'O-XYL' };

% ----------------------------------------------------------------------- %
% Database
% ----------------------------------------------------------------------- %        
filename = 'Database.xml';
foldername = '../data/';
components_database = PopulateDatabase(foldername, filename);
[database] = ReshapeDatabase(components_database, family_names, list_species);

% ----------------------------------------------------------------------- %
% [UD] Random number generator 
% ----------------------------------------------------------------------- %
rng_seed = 0;               % initial seed
rng_generator = 'twister';  % 'twister' | 'combRecursive' | 'multFibonacci'

% ----------------------------------------------------------------------- %
% [UD] Target
% ----------------------------------------------------------------------- %

target.HC  = 1.91;          % H/C ratio (experiments)
target.MW  = 153.;          % Molecular weight (kg/kmol experiments)
target.CN  = 43.;           % Cetane number (experiments)
target.TSI = 50.;           % Threshold Soot Index (experiments)
target.mu = 1;              % Viscosity (cP @298K, experiments)
target.YSI = 10;            % Yield Sooting Index (experiments)

% Viscosity
target.Tvis = 298.;         % Temperature at which viscosity is evaluated (in K)

% Distillation curve
target.vol = ReadFromFileTable('../data/DistillationCurve.JetA1', 1);
target.Td  = ReadFromFileTable('../data/DistillationCurve.JetA1', 2);
target.P = 760.;        % distillation curve: pressure (mmHg)
target.delta = 1.;      % distillation curve: user-defined step (%)

% Standard deviations (exp. estimated)
target.sigma_MW  = 10.;          % Molecular weight (in kg/kmol)
target.sigma_HC  = 0.1;          % H/C ratio 
target.sigma_CN  = 2.5;          % Cetane number
target.sigma_TSI = 15.;          % Threshold Soot Index
target.sigma_DC  = 2.5;          % Distillation temperature (in C)
target.sigma_mu = 1.;            % Viscosity (in cP)
target.sigma_YSI  = 10.;         % Yield Sooting Index (experiments)

% Additional options: genetic algorithm only
target.global_constraints_nspecies = true;      % constraints on total number of species
target.family_constraints_nspecies = true;      % constraints on number of species per family 
target.ln = [1 0 1 1];                          % minimum number of species per family
target.un = [2 0 1 1];                          % maximum number of species per family
target.nmax = 4;                                % maximum number of species (total)
target.minx = 0.;                               % minimum mole fractions for individual component
target.v_constraints = false;                   % constraints on min/max amounts of families (volume fractions)
target.lv = [0.00 0.00 0.00 0.00];              % min amounts of families (volume fractions)
target.uv = [1.00 1.00 1.00 1.00];              % max amounts of families (volume fractions)

% ----------------------------------------------------------------------- %
% [UD] Optimizer
% ----------------------------------------------------------------------- %

% Error type function
optimizer.error_type = 2;
    
% Objective function weights
optimizer.weight_MW  = 0.;          % molecular weight
optimizer.weight_HC  = 0.;          % H/C ratio
optimizer.weight_CN  = 0.;          % cetane number
optimizer.weight_TSI = 0.;          % threshold soot index
optimizer.weight_mu  = 0.;          % viscosity
optimizer.weight_YSI = 0.;          % yield soot index
optimizer.weight_DC  = 1.;          % distillation curve


% Absolute constraints (ga only)
optimizer.abs_constraints = false;  % absolute contraints
optimizer.delta_abs_MW  = 10.0;     % molecular weight (kg/kmol)
optimizer.delta_abs_HC  = 0.20;     % H/C ratio (-)
optimizer.delta_abs_CN  = 10.0;     % cetane number (0-100)
optimizer.delta_abs_TSI = 100.;     % threshold soot index (0-100)
optimizer.delta_abs_mu  = 1.;       % viscosity (cP)
optimizer.delta_abs_YSI = 100.;     % yield soot index

% Relative constraints (ga only)
optimizer.rel_constraints = false;  % relative contraints
optimizer.delta_rel_MW  = 0.025;    % molecular weight
optimizer.delta_rel_HC  = 0.025;    % H/C ratio
optimizer.delta_rel_CN  = 0.025;    % cetane number
optimizer.delta_rel_TSI = 0.025;    % threshold soot index
optimizer.delta_rel_mu  = 0.025;    % viscosity
optimizer.delta_rel_YSI = 0.025;    % yield soot index

% Additional options: genetic algorithm only
genetic_algorithm.population_size = 200;        % population size
genetic_algorithm.crossover_fraction = 0.20;    % cross-over fraction (low means more exploration)

% Additional options: particle swarm algorithm only
particle_swarm.swarm_size = 200;    % swarm size

% Additional options: pattern search algorithm only
pattern_search.poll_method = 'GPSPositiveBasis2N';  % pattern the algorithm uses to create the mesh: 'GPSPositiveBasis2N' | 'GSSPositiveBasis2N' | 'MADSPositiveBasis2N'        
pattern_search.use_complete_poll = true;            % whether all the points in the current mesh must be polled at each iteration
pattern_search.poll_order_algorithm = 'Random';     % order in which algorithm searches points in current mesh: 'Random' | 'Success' | 'Consecutive'
pattern_search.search_fct = 'searchneldermead';     % optional search step: 'None' | 'searchga' | 'searchlhs' | 'searchneldermead'


% ----------------------------------------------------------------------- %
% [UD] First guess solution
% ----------------------------------------------------------------------- %
x_fg = {    
            'NC12H22',  0.50; 
            'NC14H30',  0.20; 
            'MCH',      0.05;
            'O-XYL',    0.25;
       };

   
% ----------------------------------------------------------------------- %
% Distillation curves for all the species in the database
% ----------------------------------------------------------------------- %
plot_distillation_curves = false;
if (plot_distillation_curves == true)
    
    rhos = ComponentDensity(database.rhoType, database.rhoCoeffs);
    figure;
    plot (target.vol, target.Td, 'o');
    xlabel('volume recovery (%)');
    ylabel('temperature (C)');
    hold on;
    for i=2:database.nx
        x = zeros(1,database.nx);
        x(1) = 0.0001;
        x(i) = 0.9999;
        [vol, Td] = DistillationCurve(target.P, x, target.delta, rhos, database.MW, database.vpType, database.vpCoeffs);
        plot (vol,Td);
        hold on;
    end
    hold off;
    
end


% ----------------------------------------------------------------------- %
% Processing input
% ----------------------------------------------------------------------- %

% Interpolation points for distillation curves
target.voli = target.vol';
target.voli(1) = target.voli(1)+0.001;
target.voli(end) = target.voli(end)-0.001;

% Sets the random number generator
rng(rng_seed, rng_generator);

% Check first-guess solution
x0 = zeros(1,database.nx);
for i=1:size(x_fg, 1)
    for j=1:database.nx
        if (strcmp(x_fg{i,1}, database.name(j)))
            x0(j) = x_fg{i,2};
            break;
        end
    end
end


% ----------------------------------------------------------------------- %
% Optimization
% ----------------------------------------------------------------------- %
if (strcmp(optimizer.type,'ga'))
    
    % Optimizer variables
    optimizer.n = 2*database.nx;            % total number of variables
    optimizer.ub = ones(1,optimizer.n);     % upper boundaries
    optimizer.ub(database.nx+1:end) = Inf;  % upper boundaries
    optimizer.lb = zeros(1,optimizer.n);    % lower boundaries
    optimizer.intcon = 1:1:database.nx;     % integer unknowns
    
    % Options
    optimizer.opts = optimoptions('ga',     'MaxStallGenerations',50,...
                                            'FunctionTolerance',1e-6,...
                                            'PopulationSize', genetic_algorithm.population_size, ...
                                            'MaxGenerations',300,...
                                            'CrossoverFraction',genetic_algorithm.crossover_fraction,...
                                            'OutputFcn', @OutputFunctionGeneticAlgorithm,...
                                            'PlotFcn',{@gaplotbestf,@gaplotstopping});

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
    [sol,fval,exitflag,output] =    ga(@ObjectiveFunction,optimizer.n, A, b, [],[], ...
                                    optimizer.lb,optimizer.ub, ...
                                    @NonLinearConstraints,optimizer.intcon, ...
                                    optimizer.opts);

    % Solution
    solution.n = sol(1:database.nx).*sol((database.nx+1):2*database.nx);
    sumn = sum(solution.n);
    solution.x = solution.n/sumn;

else
    
    optimizer.n = database.nx;
    optimizer.lb = zeros(1,optimizer.n);

    if (strcmp(optimizer.type,'particleswarm'))
        
        % Options
        optimizer.opts = optimoptions(  'particleswarm', ...
                                        'SwarmSize', particle_swarm.swarm_size, ...
                                        'OutputFcn', @OutputFunctionParticleSwarm,...
                                        'PlotFcn', @pswplotbestf);
        
        % Set 20 individuals as row vectors, the rest of the swarm is random
        optimizer.opts.InitialSwarmMatrix = repmat(x0,20,1);   

        % Optimize
        [sol,fval,exitflag,output] =    particleswarm(@ObjectiveFunction, ...
                                        optimizer.n, optimizer.lb,[], ...
                                        optimizer.opts);
                                    
    elseif (strcmp(optimizer.type,'patternsearch'))
        
        % Options
        optimizer.opts = optimoptions(  'patternsearch', ...
                                        'OutputFcn', @OutputFunctionPatternSearch,...
                                        'PollMethod', pattern_search.poll_method,...
                                        'UseCompletePoll', pattern_search.use_complete_poll,...
                                        'PollOrderAlgorithm', pattern_search.poll_order_algorithm, ...  
                                        'SearchFcn', pattern_search.search_fct, ...
                                        'PlotFcn', @psplotbestf);
        
        % Optimize
        [sol,fval,exitflag,output] =    patternsearch(@ObjectiveFunction, ...
                                        x0, [], [], [], [], optimizer.lb, [], [], ...
                                        optimizer.opts);
        
    end

    % Solution
    solution.n = sol(1:database.nx);
    sumn = sum(solution.n);
    solution.x = solution.n/sumn;
    
end



% Molecular weight
solution.MW = sum(solution.x.*database.nH)*1+sum(solution.x.*database.nC)*12;

% H/C ratio
solution.HC = sum(solution.x.*database.nH)/sum(solution.x.*database.nC);

% Mass fractions
solution.omega = solution.x.*database.MW/solution.MW;

% Volumetric fractions (kg/m3)
solution.V = solution.omega./ComponentDensity(database.rhoType, database.rhoCoeffs);
sumV = sum(solution.V);
solution.V = solution.V/sumV;

% Fractions for families
solution.x_families = database.family_table*solution.x';
solution.omega_families = database.family_table*solution.omega';
solution.V_families = database.family_table*solution.V';
    
% Cetane number
solution.CN = solution.V*database.CN';

% TSI
solution.TSI = solution.x*database.TSI';

% Distillation curve
rhos = ComponentDensity(database.rhoType, database.rhoCoeffs);
[vol, Td] = DistillationCurve(target.P, solution.x, target.delta, rhos, database.MW, database.vpType, database.vpCoeffs);
plot (vol,Td, '-', target.vol, target.Td, 'o');
xlabel('volume recovery (%)');
ylabel('temperature (C)');


fprintf('------------------------------------------------------------------\n');
fprintf('   Optimal Surrogate Composition: individual species              \n');
fprintf('------------------------------------------------------------------\n');
fprintf('   Index Family          Name          x      omega          V    \n');
for i=1:database.nx
    if (solution.x(i) > 0.)
        fprintf(' * [%2d]  [%1d]  %15s %10.4f %10.4f %10.4f\n', ...
            i, database.family_index(i), database.name(i), ...
            solution.x(i)*100, solution.omega(i)*100, solution.V(i)*100);
    end
end
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

fprintf('------------------------------------------------------------------\n');
fprintf('   Optimal Surrogate Composition: families                        \n');
fprintf('------------------------------------------------------------------\n');
fprintf('   Family          Name          x      omega          V          \n');
for i=1:length(family_names)
    fprintf(' * [%1d]  %15s %10.4f %10.4f %10.4f\n', i, family_names{i}, ...
        solution.x_families(i)*100, solution.omega_families(i)*100, solution.V_families(i)*100);
end
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

fprintf('------------------------------------------------------------------\n');
fprintf('   Target values                                                  \n');
fprintf('------------------------------------------------------------------\n');
fprintf(' * HC:  Target =%8.3f - Solution =%8.3f\n', target.HC,  solution.HC);                    
fprintf(' * MW:  Target =%8.3f - Solution =%8.3f\n', target.MW,  solution.MW); 
fprintf(' * CN:  Target =%8.3f - Solution =%8.3f\n', target.CN,  solution.CN);
fprintf(' * TSI: Target =%8.3f - Solution =%8.3f\n', target.TSI, solution.TSI);
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

if (strcmp(optimizer.type,'ga'))
    fprintf('------------------------------------------------------------------\n');
    fprintf('   Optimization: summary                                          \n');
    fprintf('------------------------------------------------------------------\n');
    fprintf(' * Number of generations:          %d\n', output.generations);
    fprintf(' * Number of function evaluations: %d\n', output.funccount);
    fprintf(' * Best function value found:      %g\n', fval);
    fprintf('------------------------------------------------------------------\n');
end
