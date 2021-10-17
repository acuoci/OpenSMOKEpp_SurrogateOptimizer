function solution = runPS(allowedspecies,x0,optimizer,target,pattern_search,database,net,idtkriging,lbvmodel,isimu,totsimu)
   if nargin < 10, isimu = 0; totsimu = 0; end
    optimizer.type='patternsearch';
    optimizer.opts = optimoptions(  'patternsearch', ...
                     'PollMethod', pattern_search.poll_method,...
                     'UseCompletePoll', pattern_search.use_complete_poll,...
                     'PollOrderAlgorithm', pattern_search.poll_order_algorithm, ...  
                     'SearchFcn', pattern_search.search_fct, ...
                     'UseVectorized',false, ...
                     'UseParallel',false,...
                     'Display','off');
    optimizer.n = database.nx;
    optimizer.lb = zeros(1,nnz(allowedspecies));
    optimizer.ub = optimizer.lb + 1;
    a = find(allowedspecies);
    xFG = x0(a);
    
    % Optimize
    try
    [sol,fval,exitflag,output] = patternsearch(@(z)ObjectiveFunction(z,allowedspecies,database,target,optimizer,net,idtkriging,lbvmodel), ...
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
    solution.exitflag = NaN;
    solution.x = zeros(1,length(allowedspecies)) + NaN;
%     fprintf('Some error occured: %s\n',ex.message);
    end
    solution.solvername = 'patternsearch';
    solution.allowedspecies = allowedspecies;
    %fprintf('[patternsearch]\topt end (flag %d)\n',solution.exitflag);
    fprintf('  [%3.3d/%3.3d]\tpatternsearch end - flag %d\n',isimu,totsimu,solution.exitflag);
end