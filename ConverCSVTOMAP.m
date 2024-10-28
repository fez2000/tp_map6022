%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: Saison20232024.csv
%


%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ";";

% Specify column names and types
VariableNames = ["TempsDeJeuenSec", "ButsEsprs", "Lancers"];
 
opts.VariableNames = VariableNames;
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
TableSaison20232024 = readtable("datasets\Saison20232024.csv", opts);
Data = table2array(TableSaison20232024);

save datasets\Saison20232024.mat Data VariableNames

%% Clear temporary variables
clear opts