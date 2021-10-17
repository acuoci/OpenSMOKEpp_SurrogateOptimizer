% Optimizer for fuel surrogates > Output file

% This script is run automatically by the main script, 'OptimizeSurrogate'

fid = fopen('output.out','w');
clk = clock;
timemin = floor(timeopt/60); timesec = round(timeopt - 60*timemin);
fprintf(fid,'Output file of Surrogate Fuel optimization:\n');
fprintf(fid,' %2.2d:%2.2d of %2.2d/%2.2d/%4.4d, optimization ended after %d minutes and %d seconds\n',clk(4),clk(5),clk(3),clk(2),clk(1),timemin,timesec);

fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Optimal Surrogate Composition: individual species              \n');
fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Index Family          Name      x [%%]  omega [%%]      V [%%]\n');
for i=1:database.nx
    if (solution.x(i) > 0.)
        fprintf(fid,' * [%2d]  [%1d]  %15s %10.4f %10.4f %10.4f\n', ...
            i, database.family_index(i), database.name(i), ...
            solution.x(i)*100, solution.omega(i)*100, solution.V(i)*100);
    end
end

fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Optimal Surrogate Composition: families                        \n');
fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Family           Name     x [%%]   omega [%%]      V [%%]          \n');
for i=1:length(family_names)
    fprintf(fid,' * [%1d]  %16s %10.4f %10.4f %10.4f\n', i, family_names{i}, ...
        solution.x_families(i)*100, solution.omega_families(i)*100, solution.V_families(i)*100);
end

fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Target values                                                  \n');
fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'        Target value   Solution value    Weight    Squared relative error\n');
fprintf(fid,' * HC:  %8.3f      %8.3f           %2.1e   %5.4e\n', target.HC,  solution.HC, optimizer.weight_HC,solution.Fsolution.F_HC);  
fprintf(fid,' * MW:  %8.3f      %8.3f           %2.1e   %5.4e\n', target.MW,  solution.MW, optimizer.weight_MW,solution.Fsolution.F_MW); 
fprintf(fid,' * CN:  %8.3f      %8.3f           %2.1e   %5.4e\n', target.CN,  solution.CN, optimizer.weight_CN,solution.Fsolution.F_CN); 
fprintf(fid,' * TSI: %8.3f      %8.3f           %2.1e   %5.4e\n', target.TSI,  solution.TSI, optimizer.weight_TSI,solution.Fsolution.F_TSI); 
fprintf(fid,' * mu:  %8.3f      %8.3f           %2.1e   %5.4e\n', target.mu,  solution.mu, optimizer.weight_mu,solution.Fsolution.F_mu); 
fprintf(fid,' * YSI: %8.3f      %8.3f           %2.1e   %5.4e\n', target.YSI,  solution.YSI, optimizer.weight_YSI,solution.Fsolution.F_YSI); 
fprintf(fid,' * rho: %8.3f      %8.3f           %2.1e   %5.4e\n', target.rho,  solution.rho, optimizer.weight_rho,solution.Fsolution.F_rho);
fprintf(fid,' * DC:                                   %2.1e   %5.4e\n', optimizer.weight_DC,solution.Fsolution.F_DC);
fprintf(fid,' * IDT:                                  %2.1e   %5.4e\n', optimizer.weight_idt,solution.Fsolution.F_IDT);
fprintf(fid,' * LBV:                                  %2.1e   %5.4e\n', optimizer.weight_lbv,solution.Fsolution.F_LBV);      

fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Metamodels used                                                \n');
fprintf(fid,'------------------------------------------------------------------\n');
if target.doIDTLBV == true
fprintf(fid,' * IDT model:\t%s\n',target.IdtCalculationType);
fprintf(fid,' * LBV model:\t%s\n',lbvmodel.LbvCalculationType);
else
fprintf(fid,'   IDT and LBV are not available with this species palette\n'); 
end

fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'   Optimization summary                                           \n');
fprintf(fid,'------------------------------------------------------------------\n');
fprintf(fid,'     Try#   | Obj.Function | ExFl | Iter | Solver  | Composition found ');
fprintf(fid,' (ordered as'); for i=1:numel(database.name), fprintf(fid,' %s',database.name(i)); end, fprintf(fid,' molar fractions)\n');
    fprintf(fid,' 1>        ');
    fprintf(fid,' |              | %d    ',solution1.exitflag);
    fprintf(fid,'| %4.4d | %s      |',solution1.output.generations, solution1.solvername);
    fprintf(fid,'  %5.4f', solution1.x ); fprintf(fid,'\n');
for i = 1:length(solutions2)
    fprintf(fid,' 2> try%3.3d ',i);
    if solutions2{i}.done == true
    fprintf(fid,' | %5.5f      | %d    ',solutions2{i}.fval,solutions2{i}.exitflag);
    fprintf(fid,'| %4.4d | %s |',solutions2{i}.output.iterations, solutions2{i}.solvername);
    fprintf(fid,'  %5.4f', solutions2{i}.x );
    fprintf(fid,'\n');
    else
    fprintf(fid,' Optimization failed.\n'); 
    end
end
fprintf(fid,'   Has been chosen the try #%d\n',best);
fprintf(fid,'------------------------------------------------------------------\n');
fclose(fid);
fprintf('Results saved to ''output.out''\n');
clear clk i fid timemin timesec;