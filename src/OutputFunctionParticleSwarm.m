function stop = OutputFunctionParticleSwarm(options,state)

    global optimizer;
    global fOut;
    global database;
    
    stop = false;
    
    switch state
        case 'init'
            
            fOut = fopen('output.txt','w');
            fprintf(fOut, '%7s %12s %12s %12s %12s', '#', 'F_Obj', 'F_MW', 'F_HC', 'F_DC');
            for i=1:database.nx
                fprintf(fOut, '%12s', database.name(i));
            end
            fprintf(fOut, '\n');
            
        case 'iter'
            
            % Recover mole fractions
            x = options.bestx/sum(options.bestx);
            
            % Print on file
            fprintf(fOut, '%7d %12.3e %12.3e %12.3e %12.3e', ...
                options.iteration, options.bestfval, optimizer.F_MW, optimizer.F_HC, optimizer.F_DC);
            for i=1:database.nx
                fprintf(fOut, '%12.6f', x(i));
            end
            fprintf(fOut, '\n');
            
            % Print on the screen
            fprintf('Iteration: %d - Best objf: %f\n', options.iteration, options.bestfval);
            
        case 'done'
            
            fclose(fOut);
            
    end
end