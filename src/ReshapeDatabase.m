function database = ...
    ReshapeDatabase(components_database, family_names, list_species)

% ReshapeDatabase - Extracts from an existing database (see
% PopulateDatabase function) the properties of a list of species selcted 
% by the user which will be used in the optimization procedure 
%
% Syntax:  database = ...
%          ReshapeDatabase(components_database, family_names, list_species)
%
% Inputs:
%    components_database - vector of individual components
%    family_names - names of families
%    list_species - user defined list of species
%
% Outputs:
%    database - database of user defined selected species
%
% --------------------------- BEGIN CODE -------------------------------- %

    % Reshape database
    database.nx = length(list_species);

    % Pre-allocation
    database.family_table = zeros(length(family_names), length(list_species));
    database.family_index = zeros(1, length(list_species));

    % Loop
    count = 0;
    for j = 1:length(list_species)

        % Find species from the database
        k = 0;
        for i=1:length(components_database)
            if (strcmp(components_database(i).name, list_species{j}))
                k = i;
                break;
            end
        end
        
        % Reshape database
        if (k ~= 0)
            
            count = count + 1;
            
            database.name(count) = string(components_database(k).name);

            database.nC(count) = components_database(k).nC;
            database.nH(count) = components_database(k).nH;
            database.nO(count) = components_database(k).nO;
            database.MW(count) = components_database(k).MW;
            database.CN(count) = components_database(k).CN;
            database.TSI(count) = components_database(k).TSI;

            database.rhoType(count) = string(components_database(k).rhoType);
            database.rhoCoeffs(count) = components_database(k).rhoCoeffs;

            database.vpType(count) = string(components_database(k).distillationType);
            database.vpCoeffs{count} = components_database(k).distillationCoeffs;

            % Recognize families
            for i=1:length(family_names)
                if (strcmp(components_database(k).family,family_names{i}))
                    database.family_table(i,count) = 1;
                    database.family_index(count) = i;
                end
            end
            if (database.family_index(count) == 0)
                error('The %s family was not declared in the family_names object \n', components_database(k).family);
            end
        
        else
            error('The following component is not included in the database: %s \n', list_species{j} );
        end

    end
    
    % Normal boiling temperature (in C)
    for j = 1:length(list_species)
        x0 = zeros(1,length(list_species));
        x0(j)=1.;
        database.Tbn(j) = Flash(0, 760., x0, 0., database.vpType, database.vpCoeffs);
    end
    
    fprintf('------------------------------------------------------------------\n');
    fprintf('   Database: individual species                                   \n');
    fprintf('------------------------------------------------------------------\n');
    fprintf('       Species Family      MW     H/C     CN     TSI    Tb(C)     \n');
    for j = 1:length(list_species)
        fprintf('%14s %6d %7.2f %7.3f %6.2f  %6.2f  %7.3f\n', ...
            list_species{j}, database.family_index(j), ...
            database.MW(j), database.nH(j)/database.nC(j), database.CN(j),...
            database.TSI(j), database.Tbn(j) );
    end
    fprintf('------------------------------------------------------------------\n');