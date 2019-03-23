function components_database = ...
         PopulateDatabase(foldername, filename)
     
% PopulateDatabase - Reads a XML database containing the relevant data
% such as atomic composition, cetane number, thermodynamic properties, etc.
% for individual species.
%
% Syntax:  components_database = ...
%          PopulateDatabase(foldername, filename, components_database)
%
% Inputs:
%    foldername - Name of folder containing the database
%    filename - Name of XML file corresponding to the database
%
% Outputs:
%    components_database - vector of individual components
%
% --------------------------- BEGIN CODE -------------------------------- %

    % PARSEXML Convert XML file to a MATLAB structure.
    try
       name = strcat(foldername, filename);
       tree = xmlread(name);
    catch
       error('Failed to read XML file %s.',filename);
    end

    % Recurse over child nodes. This could run into problems with very deeply nested trees.
    try

       mainNode = tree.getDocumentElement;
       componentsNode  = mainNode.getElementsByTagName('Components');

       for k = 0:componentsNode.getLength-1

           thisComponentSet = componentsNode.item(k);   
           componentNode = thisComponentSet.getElementsByTagName('component');

           for j = 0:componentNode.getLength-1

                 thisComponent = componentNode.item(j);

                 component_name = strtrim(char(thisComponent.getAttribute('name')));

                 namesNode = thisComponent.getElementsByTagName('names').item(0);
                 names = strsplit(char(namesNode.getTextContent));

                 familyNode = thisComponent.getElementsByTagName('family').item(0);
                 family = strtrim(char(familyNode.getTextContent));

                 compositionNode = thisComponent.getElementsByTagName('composition').item(0);
                 nC = str2double(compositionNode.getAttribute('C'));
                 nH = str2double(compositionNode.getAttribute('H'));
                 nO = str2double(compositionNode.getAttribute('O'));

                 % Density
                 densityNode = thisComponent.getElementsByTagName('density').item(0);
                 densityType = strtrim(char(densityNode.getAttribute('type')));
                 if (densityType == 'constant')
                     densityCoeffs = str2double(densityNode.getTextContent);
                 else
                 end
                 
                 % Viscosity
                 viscosityNode = thisComponent.getElementsByTagName('viscosity').item(0);
                 viscosityType = strtrim(char(viscosityNode.getAttribute('type')));
                 if (viscosityType == 'law-1')
                    viscosityCoeffs = str2double(strsplit(char(viscosityNode.getTextContent)));
                    viscosityCoeffs = viscosityCoeffs(~isnan(viscosityCoeffs));
                 else
                 end
                 
                 % Vapor pressure
                 distillationNode = thisComponent.getElementsByTagName('distillation').item(0);
                 distillationType = strtrim(char(distillationNode.getAttribute('type')));
                 if (distillationType == 'Antoine3')
                     distillationCoeffs = str2double(strsplit(char(distillationNode.getTextContent)));
                     distillationCoeffs = distillationCoeffs(~isnan(distillationCoeffs));
                 else
                 end
                 
                 % Cetane number
                 cetaneNode = thisComponent.getElementsByTagName('cetane').item(0);
                 CN = str2double(strtrim(char(cetaneNode.getTextContent)));
                 
                 % Threshold soot index
                 tsiNode = thisComponent.getElementsByTagName('TSI').item(0);
                 TSI = str2double(strtrim(char(tsiNode.getTextContent)));
                 
                 % Yield soot index
                 ysiNode = thisComponent.getElementsByTagName('YSI').item(0);
                 YSI = str2double(strtrim(char(ysiNode.getTextContent)));
                 
                 
                 index = j+1;
                 components_database(index).name = component_name;
                 components_database(index).family = family;
                 components_database(index).nC = nC;
                 components_database(index).nH = nH;
                 components_database(index).nO = nO;
                 components_database(index).MW = nC*12. + nH*1. + nO*16.;
                 components_database(index).rhoType = densityType;
                 components_database(index).rhoCoeffs = densityCoeffs;
                 components_database(index).muType = viscosityType;
                 components_database(index).muCoeffs = viscosityCoeffs;
                 components_database(index).distillationType = distillationType;
                 components_database(index).distillationCoeffs = distillationCoeffs;
                 components_database(index).CN = CN;
                 components_database(index).TSI = TSI;
                 components_database(index).YSI = YSI;

             end

       end

    catch
       error('Unable to parse XML file %s.', filename);
    end
    