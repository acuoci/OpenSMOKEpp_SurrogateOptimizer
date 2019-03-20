function y = ReadFromFileTable(filename, col)

% ReadFromFileTable - Reads a single column of data from ASCII file
%
% Syntax:  ReadFromFileTable(filename, col)
%
% Inputs:
%    filename - file name (including relative or absolute path)
%    col - index of column to read
%
% Outputs:
%    y - column of data
%
% --------------------------- BEGIN CODE -------------------------------- %

    A = importdata(filename);
    y = A(:,col);
