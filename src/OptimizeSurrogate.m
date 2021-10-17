% ----------------------------------------------------------------------- %
%                   OPTIMIZER FOR FUEL SURROGATES                         %
% ----------------------------------------------------------------------- %
% Simone Pertesana, Alberto Cuoci, Marco Mehl                             %
% CRECK Modeling Lab                                                      %
% Department of Chemistry, Materials, and Chemical Engineering            %
% Politecnico di Milano (Italy)                                           %
% ----------------------------------------------------------------------- %

clear variables;
clear global;
close all;

global target;          % definition of targets for optimization
global database;        % database of individual components
global optimizer;       % optimizer options
global net              % neural network options
global idtkriging;      % kriging model (for IDT)
global lbvmodel;        % LBV polynomial model

% ----------------------------------------------------------------------- %
% Species available in the database
% ----------------------------------------------------------------------- %
% n-alkanes:    'NC4H10'	'NC5H12'	'NC7H16'	'NC10H22'	'NC12H26'	'NC16H34'	
% i-alkanes:    'IC5H12'	'IC8H18'	'IC12H26'	'IC16H34'	
% c1-alkanes:   'CHX'       'MCYC6'     'DECALIN'	
% aromatics:    'TMBENZ'	'C6H5C4H9'	'C10H7CH3'	'O-XYL'/'XYLENE'
% naphto-aromatics:         'TETRA'
% ----------------------------------------------------------------------- %
% Species available for the IDT and LBV curves
% ----------------------------------------------------------------------- %
% 'IC16H34' 'NC10H22' 'NC12H26' 'NC7H16' 'TMBENZ' 'XYLENE'

showPlots  = false; % Display a plot comparing target and surrogate: 'true' | 'false'
saveOutput = true;  % Save an 'output.out' file: 'true' | 'false'

% ----------------------------------------------------------------------- %
% Surrogate specs
% ----------------------------------------------------------------------- %    
% List of species that can be used to build the surrogate
target.list_species   = { 'IC16H34', 'NC10H22', 'NC12H26', 'NC7H16', 'TMBENZ', 'XYLENE' };
% target.list_species = {'NC7H16','NC10H22','NC12H26','NC16H34','IC8H18','IC12H26','IC16H34','TMBENZ','C6H5C4H9','TETRA','DECALIN','MCYC6','C10H7CH3'};
target.active_species = true(1,numel(target.list_species)); 
%   It's possible to fastly disable few species setting to 0/false the
%   corresponding element of 'target.active_species'

target.atLeastOneSpeciesPerFamily = false;
target.nSpecies = 0; % Impose # of species of the surrogate - 0=unset

%  . . . . . . . . . . . . . . . . . Properties to be optimized (0=no/1=yes):
                                     % MW  H/C CN  TSI mu  YSI DC  rho IDT LBV
target.optimizedProperties = logical([ 0   1   1   0   0   0   1   1   0   0  ]);

% ----------------------------------------------------------------------- %
% Optimizer specs
% ----------------------------------------------------------------------- %    
optimizer.parallel_computing = true; % Try to explot parallel computing: 'true' | 'false'
optimizer.solver1step = 'ga';       % 'ga'
optimizer.solver2step = 'fmincon';  % 'fmincon' | 'patternsearch' | 'both'
max_species_change = 4; % # of species that can be added/removed from step 1 to step 2 in the optimization process

% ----------------------------------------------------------------------- %
% Database
% ----------------------------------------------------------------------- %        
filename = 'Database.xml';          foldername = '../data/';
[components_database,family_names] = PopulateDatabase(foldername, filename);
[database] = ReshapeDatabase(components_database, family_names, target.list_species);
clear filename foldername;

% ----------------------------------------------------------------------- %
% [UD] Random number generator 
% ----------------------------------------------------------------------- %
rng_seed = 0;               % initial seed
rng_generator = 'twister';  % 'twister' | 'combRecursive' | 'multFibonacci'

% ----------------------------------------------------------------------- %
% Neural Network (for IDTs - old code)
% ----------------------------------------------------------------------- % 
net.species = [ "DECALIN", "NC10H22", "NC12H26", "IC16H34", "TMBENZ", "XYLENE" ];
ANN=load('../data/meta-models/ANN/NetDECALIN.mat', 'net');
net.net = ANN.net;
net.logarithm = "true";     % "true" | "false"
clear ANN;

% ----------------------------------------------------------------------- %
% Kriging Model (for IDTs)
% ----------------------------------------------------------------------- % 
idtkriging.species = ["IC16H34","NC10H22","NC12H26","NC7H16","TMBENZ","XYLENE"];
    % Check if IDT and LBV are available with the current species subset
    target.doIDTLBV = true;
    for j = 1:database.nx
        specieok = false;
        if target.active_species(j) == false, specieok = true;
        else
        for i = 1:numel(idtkriging.species)
            if strcmp(idtkriging.species(i), target.list_species{j}), specieok = true; end
        end
        end
        if specieok == false, target.doIDTLBV = false; end
    end
    if target.doIDTLBV == false, warning('IDT and LBV are not available with the current species palette'); end
    clear i j specieok
if target.doIDTLBV == true
idtm = load('../data/meta-models/Kriging/IDT_Kriging_6.mat','idtmodel');
idtkriging.model = idtm.idtmodel;
idtkriging.logarithm = "true";
fprintf('IDT metamodel loaded.\n');
clear idtm;
end

% ----------------------------------------------------------------------- %
% LBV model
% ----------------------------------------------------------------------- % 
if target.doIDTLBV == true
lbvmodel.LbvCalculationType = "linear";   % 'linear' | 'quadratic' | 'cubic' | 'quadratic-mc' (quadratic with mixed coefficients)
lbvpuresfolder = '../data/pure-components/LBV_400K_1bar/';
lbvmodel.sLpures = LBVPreProcessing(target.list_species,lbvpuresfolder);
v_pures=zeros(7,length(target.list_species));  % These vectors have as first size the
phi_pures=zeros(7,length(target.list_species));% number of points of the LBV curve
for i=1:length(target.list_species)
A1=importdata([lbvpuresfolder target.list_species{i} '.out']);
v_pures(:,i) = A1.data(:,2);
phi_pures(:,i) = A1.data(:,3);
end
[~,I]=sort(phi_pures);
for j = 1:size(I,2), v_pures(:,j)=v_pures(I(:,j),j); end
lbvmodel.sLpures = v_pures;
lbvmodel.phipures = phi_pures(:,1);
if strcmp(lbvmodel.LbvCalculationType,"linear")
    lbm = load('../data/meta-models/Polynomials/LBV_1.mat');
elseif strcmp(lbvmodel.LbvCalculationType,"quadratic")
    lbm = load('../data/meta-models/Polynomials/LBV_2.mat');
elseif strcmp(lbvmodel.LbvCalculationType,"cubic")
    lbm = load('../data/meta-models/Polynomials/LBV_3.mat');
elseif strcmp(lbvmodel.LbvCalculationType,"quadratic-mc")
    lbm = load('../data/meta-models/Polynomials/LBV_2-mc.mat');
end
lbvmodel.LBVcoeffs = lbm.parameters;
lbvmodel.LBVexponents = lbm.MatCoeff;
fprintf('LBV metamodel loaded.\n');
clear lbm A1 v_pures phi_pures i I j;
else, lbvmodel = []; 
end

% ----------------------------------------------------------------------- %
% [UD] Target
% ----------------------------------------------------------------------- %
fprintf('Setting target\n');
%                                                                       JetA1 | POSF-4658 |
target.HC  = 1.957;         % H/C ratio (experiments)                   1.94  | 1.957     |
target.MW  = 142.;          % Molecular weight (kg/kmol experiments)    154   | 142       |
target.CN  = 47.1;          % Cetane number (experiments)               47    | 47.1      |
target.TSI = 50.;           % Threshold Soot Index (experiments)        50.
target.mu = 2.3878;         % Viscosity (cP @298K, experiments)    2.6672@298K| 2.3878@273K
target.YSI = 75.5;          % Yield Sooting Index (experiments)         60.   | 75.5      |
target.rho = 799.38;        % Density                                   801   | 799.38    |

% Viscosity
target.Tvis = 273.;         % Temperature at which viscosity is evaluated (in K)

% Distillation curve
dcfile = '../data/exp/DistillationCurve.POSF4658';
DC = importdata(dcfile);
target.vol = DC(:,1); 
target.Td  = DC(:,2);
target.P = 760.;         % distillation curve: pressure (mmHg)
target.delta = 10.;      % distillation curve: user-defined step (%)
fprintf('  DC  imported from ''%s''\n',dcfile);
clear DC dcfile;

% Ignition delay times
target.IdtCalculationType = 'kriging';  % 'net' | 'direct' | 'kriging'
target.tauType = 'Tslope';              % 'Tslope' | [...]
idtfile = '../data/exp/IDT.POSF4658';
A1=importdata(idtfile); 
A1.data = sortrows(A1.data,2); % Order in ascending T
target.TIdt = A1.data(:,2);             % Temperature [K];
target.Idt = A1.data(:,4);              % Ignition delay time [s]
target.PIdt = A1.data(:,3).*1.01325;    % Pressure [atm]
target.phiIdt = A1.data(:,1);           % Equivalence ratio
fprintf('  IDT imported from ''%s''\n',idtfile);
clear A1 idtfile;

% Laminar burning velocities
lbvfile = '../data/exp/LBV.POSF4658';
A2=importdata(lbvfile);
targetphiLBV = A2.data(:,1);
targetLBV = A2.data(:,2);
[target.phiLBV,I2]=sort(targetphiLBV);
for j = 1:length(I2), target.lbv(j)=targetLBV(I2(j)); end
fprintf('  LBV imported from ''%s''\n',lbvfile);
clear A2 targetphiLBV targetLBV j I2 lbvfile;

% ----------------------------------------------------------------------- %
% [UD] Optimizer options
% ----------------------------------------------------------------------- %
% Note that these constraints work only at the first level that uses 'ga'
%  and may not work with parallel computing enabled

% Additional options: genetic algorithm only
target.global_constraints_nspecies = false;   % constraints on total number of species
target.family_constraints_nspecies = false;   % constraints on number of species per family 
target.ln = [0 0 0 0];                        % minimum number of species per family
target.un = [4 4 4 4];                        % maximum number of species per family
target.nmax = 11;                             % maximum number of species (total)
target.minx = 0.;                             % minimum mole fractions for individual component
target.v_constraints = false;                 % constraints on min/max amounts of families (volume fractions) 
target.lv = [0.1  0.2  0.0  0.1];             % min amounts of families (volume fractions)   
target.uv = [1.0  1.0  1.0  1.0];             % max amounts of families (volume fractions)

% Absolute constraints (ga only)
optimizer.abs_constraints = false;  % absolute contraints
optimizer.delta_abs_MW  = 10.0;     % molecular weight (kg/kmol)
optimizer.delta_abs_HC  = 0.20;     % H/C ratio (-)
optimizer.delta_abs_CN  = 10.0;     % cetane number (0-100)
optimizer.delta_abs_TSI = 100.;     % threshold soot index (0-100)
optimizer.delta_abs_mu  = 1.;       % viscosity (cP)
optimizer.delta_abs_YSI = 100.;     % yield soot index
optimizer.delta_abs_rho = 10.;      % density
optimizer.delta_abs_idt = 1.;       % idt
optimizer.delta_abs_lbv = 1.;       % laminar burning velocities

% Relative constraints (ga only)
optimizer.rel_constraints = false;  % relative contraints
optimizer.delta_rel_MW  = 0.025;    % molecular weight
optimizer.delta_rel_HC  = 0.025;    % H/C ratio
optimizer.delta_rel_CN  = 0.025;    % cetane number
optimizer.delta_rel_TSI = 0.025;    % threshold soot index
optimizer.delta_rel_mu  = 0.025;    % viscosity
optimizer.delta_rel_YSI = 0.025;    % yield soot index
optimizer.delta_rel_rho = 0.025;    % density
optimizer.delta_rel_idt = 0.025;    % idt
optimizer.delta_rel_lbv = 0.025;    % lbv

% Additional options: genetic algorithm only
genetic_algorithm.population_size = 200;        % population size
genetic_algorithm.crossover_fraction = 0.20;    % cross-over fraction (low means more exploration)
genetic_algorithm.nparallel = 2;

% Additional options: particle swarm algorithm only
particle_swarm.swarm_size = 200;    % swarm size

% Additional options: pattern search algorithm only
pattern_search.poll_method = 'GPSPositiveBasis2N';  % pattern the algorithm uses to create the mesh: 'GPSPositiveBasis2N' | 'GSSPositiveBasis2N' | 'MADSPositiveBasis2N'        
pattern_search.use_complete_poll = true;            % whether all the points in the current mesh must be polled at each iteration
pattern_search.poll_order_algorithm = 'Random';     % order in which algorithm searches points in current mesh: 'Random' | 'Success' | 'Consecutive'
pattern_search.search_fct = 'searchneldermead';     % optional search step: 'None' | 'searchga' | 'searchlhs' | 'searchneldermead'

% ----------------------------------------------------------------------- %
% Processing input
% ----------------------------------------------------------------------- %

% Interpolation points for distillation curves
target.voli = target.vol';
target.voli(1) = target.voli(1)+0.001;
target.voli(end) = target.voli(end)-0.001;

% Sets the random number generator
rng(rng_seed, rng_generator);

% Prepare parallel pool 
if optimizer.parallel_computing == true
    try
        ParPool = gcp; % Get Current Parallel pool (if already existing - 
                       % Else, try to start a new one)
        if isempty(ParPool)
            optimizer.parallel_computing = false;
        else
            fprintf('Parallel computing will be exploited.\n');
        end
    catch
        fprintf('Parallel computing won''t be exploited - some problems were found with parallel pool\n');
        optimizer.parallel_computing = false;
    end
else
    fprintf('Parallel computing won''t be exploited\n');
end
tic;

% ----------------------------------------------------------------------- %
% Optimization #1 with ga
% ----------------------------------------------------------------------- %
optimizer.type = 'ga';
% Objective function weights
optimizer.weight_MW  = 1.;          % molecular weight
optimizer.weight_HC  = 1.;          % H/C ratio
optimizer.weight_CN  = 1.;          % cetane number
optimizer.weight_TSI = 0.;          % threshold soot index
optimizer.weight_mu  = 1.;          % viscosity
optimizer.weight_YSI = 1.;          % yield soot index
optimizer.weight_DC  = 0.;          % distillation curve
optimizer.weight_rho = 1.;          % density
optimizer.weight_idt = 0.;          % ignition curve
optimizer.weight_lbv = 0.;          % lamninar burning velocities
CheckWeights(target.optimizedProperties,target.doIDTLBV);
    
fprintf(' > Optimization step #1 (ga) ...\n');
% Loose specs
genetic_algorithm.maxStallGen = 15;
genetic_algorithm.maxGen = 100;  % 100
genetic_algorithm.FunTol = 1e-6; % 1e-6
% Run GA optimizer
solution1 = RunGA(optimizer,target,genetic_algorithm,database);
fprintf('\tfval = %e\n\tflag = %d\n\t# species = %d\n',solution1.fval,...
        solution1.exitflag,nnz(solution1.x));

% fix the species
allowedspecies = (solution1.x ~= 0);
x0 = solution1.x;
[~,solution1.Fsolution] = ...
    ObjectiveFunction([solution1.x, ones(1,length(solution1.x))],...
    true(1,length(solution1.x)),database,target,optimizer,net,idtkriging,lbvmodel);

% ----------------------------------------------------------------------- %
% Optimization #2 with patternsearch / fmincon
% ----------------------------------------------------------------------- %
% Objective function weights
optimizer.weight_MW  = 5.;         % molecular weight
optimizer.weight_HC  = 10;          % H/C ratio
optimizer.weight_CN  = 10;          % cetane number
optimizer.weight_TSI = 0.;          % threshold soot index
optimizer.weight_mu  = 5.;          % viscosity
optimizer.weight_YSI = 5.;          % yield soot index
optimizer.weight_DC  = 1.;          % distillation curve
optimizer.weight_rho = 10;          % density
optimizer.weight_idt = 0.;%1.5;     % ignition curve
optimizer.weight_lbv = 0.;%0.3      % lamninar burning velocities

optimizer.type='patternsearch';
optimizer.n = database.nx;
CheckWeights(target.optimizedProperties,target.doIDTLBV)

AllSpeciesMatrix = DoSpeciesCombinations(database,...
    target.atLeastOneSpeciesPerFamily,allowedspecies,max_species_change,target.nSpecies);
    % Starting from species found by ga, add/remove species following
    % max_species_change constraint
if isempty(AllSpeciesMatrix), AllSpeciesMatrix = allowedspecies; end % If DoSpeciesCombinations fails, use just the species found by ga

AllSpeciesMatrix = unique(and(AllSpeciesMatrix,target.active_species),'rows');
    % Then, keep only the combinations that respect the user's species
    
Nsimulations2lvl = size(AllSpeciesMatrix,1);
if strcmp(optimizer.solver2step,'both')
    Nsimulations2lvl = 2*Nsimulations2lvl;     
    AllSpeciesMatrix = repmat(AllSpeciesMatrix,2,1);
end
solutions2 = cell(1,Nsimulations2lvl);

% Make local copies of global variables, to use 'parfor'
lopt=optimizer;ltgt=target;ldb=database;lnet=net;lidt=idtkriging;llbv=lbvmodel;
secondsolver = optimizer.solver2step;

if strcmp(secondsolver,'both')
    fprintf(' > Optimization step #2 (patternsearch/fmincon, %d trials) ...\n',Nsimulations2lvl);
elseif strcmp(secondsolver,'patternsearch')
    fprintf(' > Optimization step #2 (patternsearch, %d trials) ...\n',Nsimulations2lvl);
elseif strcmp(secondsolver,'fmincon')
    fprintf(' > Optimization step #2 (fmincon, %d trials) ...\n',Nsimulations2lvl);
end

if optimizer.parallel_computing == true
% Run simulations in parallel
parfor ii = 1:Nsimulations2lvl
    if strcmp(secondsolver,'both')
        if ii <= Nsimulations2lvl/2
          solutions2{ii}=RunPS (AllSpeciesMatrix(ii,:),x0,lopt,ltgt,pattern_search,ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
        else
          solutions2{ii}=RunFMC(AllSpeciesMatrix(ii,:),x0,lopt,ltgt,[],            ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
        end
    elseif strcmp(secondsolver,'patternsearch')
          solutions2{ii}=RunPS (AllSpeciesMatrix(ii,:),x0,lopt,ltgt,pattern_search,ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
    elseif strcmp(secondsolver,'fmincon')
          solutions2{ii}=RunFMC(AllSpeciesMatrix(ii,:),x0,lopt,ltgt,[],            ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
    else
        error("'optimizer.solver2step' must be 'patternsearch', 'fmincon', or 'both'");
    end
end
else
% Else, run them in series
for ii = 1:Nsimulations2lvl
    if strcmp(secondsolver,'both')
        if ii <= Nsimulations2lvl/2
          solutions2{ii}=RunPS (AllSpeciesMatrix(ii,:),x0,lopt,ltgt,pattern_search,ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
        else
          solutions2{ii}=RunFMC(AllSpeciesMatrix(ii,:),x0,lopt,ltgt,[],            ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
        end
    elseif strcmp(secondsolver,'patternsearch')
          solutions2{ii}=RunPS (AllSpeciesMatrix(ii,:),x0,lopt,ltgt,pattern_search,ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
    elseif strcmp(secondsolver,'fmincon')
          solutions2{ii}=RunFMC(AllSpeciesMatrix(ii,:),x0,lopt,ltgt,[],            ldb,lnet,lidt,llbv,ii,Nsimulations2lvl);
    else
        error("'optimizer.solver2step' must be 'patternsearch', 'fmincon', or 'both'");
    end
end
end
clear lopt ltgt ldb lnet lidt llbv

fvals = ones(Nsimulations2lvl,1).*NaN;
for i = 1:Nsimulations2lvl
    if solutions2{i}.done == false
        solutions2{i}.fval = Inf; 
    end
    fvals(i) = solutions2{i}.fval;
end
% Take as surrogate the try with the min objective function
best = find(fvals == min(fvals)); best = best(1);
solution = solutions2{best};
fprintf('Optimal fval: %f [try %2.2d]\n',solution.fval,best);
timeopt = toc;

[~,solution.Fsolution] = ...
    ObjectiveFunction(solution.x,true(1,length(solution.allowedspecies)),...
                      database,target,optimizer,net,idtkriging,lbvmodel);
clear fvals i;

% ----------------------------------------------------------------------- %
% Post-processing
% ----------------------------------------------------------------------- %
% Compute properties of the surrogate, do plots

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

% viscosity
mus = ComponentViscosity(target.Tvis, database.muType, database.muCoeffs); 
solution.mu = sum((solution.omega).*(mus.^(1/3)))^3; % Blending rule: mu^1/3 = sum (omega_i*mu_i^1/3)

% YSI
solution.YSI = solution.omega*database.YSI';

% Density
rhos = ComponentDensity(database.rhoType, database.rhoCoeffs);
solution.rho=solution.x*rhos';

% Distillation curve
rhos = ComponentDensity(database.rhoType, database.rhoCoeffs);
[vol, Td] = DistillationCurve(target.P, solution.x, target.delta, rhos, database.MW, database.vpType, database.vpCoeffs);
solution.vol = vol; solution.Td = Td;

if target.doIDTLBV == true 
% Ignition curve
if target.IdtCalculationType == "kriging"
solution.idt = IgnitionDelayTimes(target.list_species,target.TIdt, target.PIdt, target.phiIdt, solution.x,...
   idtkriging.species, idtkriging.model, target.tauType, target.IdtCalculationType, idtkriging.logarithm); 
else
solution.idt = IgnitionDelayTimes(target.list_species, target.TIdt, target.PIdt, target.phiIdt, solution.x,...
    net.species, net.net, target.tauType, target.IdtCalculationType, net.logarithm); 
end
end
if target.doIDTLBV == true 
% laminar burning velocities
        LBVpures = lbvmodel.sLpures;
        LBVtarget = target.lbv;
        phis = target.phiLBV;
        if numel(lbvmodel.phipures) ~= numel(target.phiLBV)
            % if phi do not coincide, do an interpolation
            phimin = max([min(lbvmodel.phipures),min(target.phiLBV)]);
            phimax = min([max(lbvmodel.phipures),max(target.phiLBV)]);
            idxphi = find( (lbvmodel.phipures >= phimin) .* (lbvmodel.phipures <= phimax) );
            phis = lbvmodel.phipures(idxphi);
            LBVpures = lbvmodel.sLpures(idxphi,:);
            LBVtarget = interp1(target.phiLBV,target.lbv,phis);
        end
solution.phiLBV = phis;
solution.lbv = LaminarBurningVelocitiesPoly(LBVpures, ...
               lbvmodel.LBVcoeffs, lbvmodel.LBVexponents, solution.x);

end
% Objective function composition:
dataF=[ optimizer.weight_HC*solution.Fsolution.F_HC   , ...
        optimizer.weight_MW*solution.Fsolution.F_MW   , ...
        optimizer.weight_CN*solution.Fsolution.F_CN   , ...
        optimizer.weight_TSI*solution.Fsolution.F_TSI , ...
        optimizer.weight_mu*solution.Fsolution.F_mu   , ...
        optimizer.weight_YSI*solution.Fsolution.F_YSI , ...        
        optimizer.weight_DC*solution.Fsolution.F_DC   , ...
        optimizer.weight_rho*solution.Fsolution.F_rho , ...
        optimizer.weight_idt*solution.Fsolution.F_IDT , ...
        optimizer.weight_lbv*solution.Fsolution.F_LBV ];
dataF = dataF./sum(dataF);

if showPlots == true
    drawPlots(target,database,solution)
end

% ----------------------------------------------------------------------- %
% Display Results
% ----------------------------------------------------------------------- %
fprintf('------------------------------------------------------------------\n');
fprintf('   Target values                                                  \n');
fprintf('        Target value   Solution value   Weight    Obj.Func. contribution\n');
fprintf(' * HC:  %8.3f      %8.3f           %4.2f     %5.4e\n', target.HC,  solution.HC, optimizer.weight_HC,solution.Fsolution.F_HC);  
fprintf(' * MW:  %8.3f      %8.3f           %4.2f     %5.4e\n', target.MW,  solution.MW, optimizer.weight_MW,solution.Fsolution.F_MW); 
fprintf(' * CN:  %8.3f      %8.3f           %4.2f     %5.4e\n', target.CN,  solution.CN, optimizer.weight_CN,solution.Fsolution.F_CN); 
fprintf(' * TSI: %8.3f      %8.3f           %4.2f     %5.4e\n', target.TSI,  solution.TSI, optimizer.weight_TSI,solution.Fsolution.F_TSI); 
fprintf(' * mu:  %8.3f      %8.3f           %4.2f     %5.4e\n', target.mu,  solution.mu, optimizer.weight_mu,solution.Fsolution.F_mu); 
fprintf(' * YSI: %8.3f      %8.3f           %4.2f     %5.4e\n', target.YSI,  solution.YSI, optimizer.weight_YSI,solution.Fsolution.F_YSI); 
fprintf(' * rho: %8.3f      %8.3f           %4.2f     %5.4e\n', target.rho,  solution.rho, optimizer.weight_rho,solution.Fsolution.F_rho);
fprintf(' * DC:                                   %4.2f     %5.4e\n', optimizer.weight_DC,solution.Fsolution.F_DC);
fprintf(' * IDT:                                  %4.2f     %5.4e\n', optimizer.weight_idt,solution.Fsolution.F_IDT);
fprintf(' * LBV:                                  %4.2f     %5.4e\n', optimizer.weight_lbv,solution.Fsolution.F_LBV);    
fprintf('------------------------------------------------------------------\n');
fprintf('   Optimal Surrogate Composition: individual species              \n');
fprintf('   Index Family          Name          x      omega          V    \n');
for i=1:database.nx
    if (solution.x(i) > 0.)
        fprintf(' * [%2d]  [%1d]  %15s %10.4f %10.4f %10.4f\n', ...
            i, database.family_index(i), database.name(i), ...
            solution.x(i)*100, solution.omega(i)*100, solution.V(i)*100);
    end
end

fprintf('------------------------------------------------------------------\n');
fprintf('   Optimal Surrogate Composition: families                        \n');
fprintf('   Family           Name          x     \n');
for i=1:length(family_names)
  if solution.x_families(i) > 0
    fprintf(' * [%1d]  %16s %10.4f \n', i, family_names{i}, solution.x_families(i)*100);
  end
end
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

% Export to 'output.out'
if saveOutput == true
    WriteOutput;
end
