function y = ReadFromFileTable(filename, col)

    A = importdata(filename);
    y = A(:,col);
