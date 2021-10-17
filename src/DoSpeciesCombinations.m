function MatAllowSp = DoSpeciesCombinations(database,onePerFamily,ASin,maxchange,nspecies)
    % Returns a matrix of allowable species, changing the active species
    % of a family of species by maxchange
    
%     % if AllowedSpeciesIn respects "at least one species per family", add
%     % this constraint
%     onePerFamily = true;
%     for i = 1:size(database.family_table,1)
%         if ~any(ASin(database.family_table(i,:)==1)), onePerFamily = false; end
%     end
    
    NspPerFamily = sum(database.family_table,2);
    for i = 1:size(database.family_table,1)
    IDXfamily{i} = find(database.family_table(i,:));
    end
    
    repeat = 1;
    for currentMaxChange = min(maxchange,numel(ASin)):numel(ASin)
    repeat = repeat-1;
    
    V = [false,true];
    A = combvec(V,V);
    for i = 1:(database.nx-2)
        A = combvec(A,V);
    end
    A = A';
    ok = true(1,size(A,1));
    
    for i = 1:length(ok)
       difference = (A(i,:) ~= ASin); 
       if nnz(difference) > currentMaxChange
           ok(i) = false;
       end
       if nspecies ~=0 && nnz(A(i,:))~= nspecies
           ok(i) = false;
       end
       for ifam = 1:size(database.family_table,1) % At least one specie per family
           b = A(i,:);
           if nnz(b(IDXfamily{ifam})) == 0 && any(database.family_table(ifam,:)) && onePerFamily==true
               ok(i) = false;
           end
       end
    end
    MatAllowSp = A(ok,:);
    
    if numel(MatAllowSp) > 0
        if currentMaxChange ~= maxchange
            warning('I was not able to find a set of species that satisfies your number of species constraint. I will try changing more species than what you asked (%d instead of %d)',currentMaxChange,maxchange)
        end
        break;
    end
    end
end