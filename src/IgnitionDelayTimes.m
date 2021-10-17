function idt = ...
    IgnitionDelayTimes(species, T, P, phi, x0, net_species, net, tauType, calculationType, logarithm)
% IgnitionDelayTimes - Calculates the ignition times curve
%
% Syntax:  idt = ...
%          IgnitionDelayTimes(T, P, phi, x0, calculationType)
%
% Inputs:
%    T - temperature (in K)
%    P - pressure (in atm)
%    phi - equivalence ratio
%    x0 - mixture composition (in molar fractions)
%    calculationType - type of calculation: use 'direct' or 'net' or
%                      'kriging'
%   
% Output:
%    idt - vector of ignition delay times (ms)
%
% If you choose 'net', pass as 'net' the Artificial Neural Network, while
% if you choose 'kriging', pass as 'net' the Kriging model struct
% --------------------------- BEGIN CODE -------------------------------- %

if calculationType == "direct"
    
    list_species = string(species);
    
    line_pressure = 35;
line_equivalence_ratio = 36;
line_fuel_moles = 37;

line_temperatures = 53;

input_file = regexp( fileread('../data/idt.template.dic'), '\n', 'split');

string_pressures = strcat('\t@Pressure\t\t', " ", num2str(P), "  atm;");

string_eqRatios = strcat('\t@EquivalenceRatio\t\t', " ", num2str(phi), ";");

string_temperatures = '\t@ListOfValues\t\t';
for i=1:length(T)
    string_temperatures = strcat(string_temperatures, " ", num2str(T(i)));
end
string_temperatures = strcat(string_temperatures, " K;");

string_fuel_moles = strcat(list_species(1), " ", string(x0(1)));
for i=2:length(list_species)
    string_fuel_moles = strcat(string_fuel_moles, " ", list_species(i), " ", string(x0(i)));
end
string_fuel_moles_label = '\t@FuelMoles\t\t\t';
string_fuel_moles_label = strcat(string_fuel_moles_label, string_fuel_moles, " ;");

input_file{line_pressure} = sprintf(string_pressures);
input_file{line_equivalence_ratio} = sprintf(string_eqRatios);
input_file{line_fuel_moles} = sprintf(string_fuel_moles_label);
input_file{line_temperatures} = sprintf(string_temperatures);

fid = fopen('input.dic', 'w');
fprintf(fid, '%s\n', input_file{:});
fclose(fid);
    
    [status, output] = system('\Users\ddeme\OpenSMOKE++Suite\bin\OpenSMOKEpp_BatchReactor.exe');
    %system('"\Users\ddeme\OpenSMOKE++Suite\bin\OpenSMOKEpp_BatchReactor.exe" --input "idtpp.dic" ');

    keySet= {'Tslope', 'Pslope', 'OHslope', 'OHmax', 'CHslope', 'CHmax'};
    valueSet = [10 11 12 13 14 15];
    M = containers.Map(keySet,valueSet);

    data = importdata('IDT/ParametricAnalysisIDT.out');
    idt = data.data(:,M(tauType));
    
elseif calculationType == "net"
    
    net_species = sort(net_species);
    
    % find indeces for which strings are different
    idx = [];
    for i=1:length(species)
        a = species(i)==net_species;
        idx = [idx find(a==1)];
    end
    
    x = zeros(1, length(net_species));
    x(idx) = x0;
    
    taus_net = zeros(1,length(T));
    for i=1:length(T)
        taus_net(i) = ( net( [x T(i)]' ) );
    end
    
    if logarithm == "true"
        idt = exp(taus_net);
    elseif logarithm == "false"
        idt = taus_net;
    else
        error('Please specify a valid value for net.logarithm ("true" or "false").')
    end

elseif calculationType == "kriging"    
    
    net_species = sort(net_species);
    
    % find indexes for which strings are different
    idx = [];
    for i=1:length(species)
        a = species(i)==net_species;
        idx = [idx find(a==1)];
    end
    
    x = zeros(1, length(net_species));
    x(idx) = x0;
    x = ones(length(T),1) * x; % create a matrix (1row=1simulation)
    
    taus_krig = zeros(1,length(T));
    
    taus_krig(:) = IDTKriging(T,P.*101325,phi,x,net);
    
    if logarithm == "true"
        idt = taus_krig;
    elseif logarithm == "false"
        idt = log(taus_krig); % ?? to check
    else
        error('Please specify a valid value for net.logarithm ("true" or "false").')
    end
    
else
    
    error('Please specify a valid value ("direct", "net" or "kriging") for the IDT calculation type.');
    
end

end
