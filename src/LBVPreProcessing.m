function v_pures = LBVPreProcessing(list_species,folder)
% LBVPreProcessing - Import the laminar flame speeds for the pure 
%                    components 
%
% Syntax:  v_pures = LBVPreProcessing(list_species,folder)
%
% Inputs:
%    list_species - active species in the optimization
%    folder - char vector representing the directory where SPECIESNAME.out
%             files containing LBV data are located
%
% Output:
%    v_pures - vector of laminar burning velocities of pure components (cm/s)
%
% --------------------------- BEGIN CODE -------------------------------- %

for i=1:length(list_species)
A1=importdata([folder list_species{i} '.out']);
    if i==1
        v_pures=zeros(size(A1.data,1),length(list_species));
        phi_pures=zeros(size(A1.data,1),length(list_species));
    end
v_pures(:,i) = A1.data(:,2);
phi_pures(:,i) = A1.data(:,3);
end
[~,I]=sort(phi_pures);
for j = 1:size(I,2), v_pures(:,j)=v_pures(I(:,j),j); end

end