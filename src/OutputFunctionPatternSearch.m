function [stop,options,optchanged] = OutputFunctionPatternSearch(optimvalues,options,state)

    global optimizer;
    global fOut;
    global database;
    
    stop = false;
    optchanged = false;
    
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
            x = optimvalues.x/sum(optimvalues.x);
            
            % Print on file
            fprintf(fOut, '%7d %12.3e %12.3e %12.3e %12.3e', ...
                optimvalues.iteration, optimvalues.fval, optimizer.F_MW, optimizer.F_HC, optimizer.F_DC);
            for i=1:database.nx
                fprintf(fOut, '%12.6f', x(i));
            end
            fprintf(fOut, '\n');
            
            % Print on the screen
            fprintf('Iteration: %d - Current objf: %f\n', optimvalues.iteration, optimvalues.fval);
            
        case 'done'
            
            fclose(fOut);
            
    end
end