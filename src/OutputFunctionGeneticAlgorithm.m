function [state,options,optchanged] = OutputFunctionGeneticAlgorithm(options,state,flag)

% OutputFunctionGeneticAlgorithm - Output function for genetic optimizer
%
% Syntax:  [stop,options,optchanged] = ...
%          OutputFunctionPatternSearch(optimvalues,options,state)
%
% Inputs:
%    options - options
%    state - structure containing information about the current generation
%    flag - current status of the algorithm
%
% Outputs:
%    state - structure containing information about the current generation
%    options - options as modified by the output function
%    optchanged - boolean flag indicating changes to options. To change 
%                 options for subsequent iterations, set optchanged to true
%
% --------------------------- BEGIN CODE -------------------------------- %

    global optimizer;
    global fOut;
    global database;
    
    optchanged = false;

    switch flag
        case 'init'

            fOut = fopen('optimization.out','w');
            fprintf(fOut, '%7s %12s %12s %12s %12s', '#', 'F_Obj', 'F_MW', 'F_HC', 'F_DC');
            for i=1:database.nx
                fprintf(fOut, '%12s', database.name(i));
            end
            fprintf(fOut, '\n');

        case 'iter'

            % Find the best objective function
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);
            bestf = ObjectiveFunction(bestx);
            
            % Recover mole fractions
            x = bestx(1:database.nx).*bestx(database.nx+1:end);
            x = x/sum(x);
            
            % Print on file
            fprintf(fOut, '%7d %12.3e %12.3e %12.3e %12.3e', ...
                state.Generation, bestf, optimizer.F_MW, optimizer.F_HC, optimizer.F_DC);
            for i=1:database.nx
                fprintf(fOut, '%12.6f', x(i));
            end
            fprintf(fOut, '\n');
            
            % Stop if it is low
            if (bestf <= 1e-6)
                state.StopFlag = 'y';
                fprintf('Objective function below threshold: %f\n', bestf)
            end
            
            % Print on the screen
            fprintf('Generation: %d - Best objf: %f\n', state.Generation, bestf);
            
            % Update the fraction of mutation and crossover after 25 generations.
            if state.Generation == 50
                options.CrossoverFraction = 0.8;
                optchanged = true;
            end
            
        case 'done'
            
            fclose(fOut);
    end
        
end

function [state,options,optchanged] = OutputFunctionParticleSwarm(options,state,flag)

    global optimizer;
    global fOut;
    
    optchanged = false;
    
    switch state
        case 'init'
            
            fOut = fopen('output.txt','w');
            fprintf(fOut, '%7s %12s %12s %12s %12s \n', '#', 'F_Obj', 'F_MW', 'F_HC', 'F_DC');
            
        case 'iter'
            
            % Find the best objective function
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);
            bestf = ObjectiveFunction(bestx);
            
            % Print on file
            fprintf(fOut, '%7d %12.3e %12.3e %12.3e %12.3e \n', ...
                options.iteration, bestf, optimizer.F_MW, optimizer.F_HC, optimizer.F_DC);
            
            % Print on the screen
            fprintf('Generation: %d - Best objf: %f\n', options.iteration, bestf);
            
        case 'done'
            
            fclose(fOut);
            
    end
end