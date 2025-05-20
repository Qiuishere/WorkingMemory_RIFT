function tabdata    = make_table(matrix)
% function for creating tables based on matlab version
if nargin==0
    if verLessThan('matlab', '8.3')
    tabdata   = dataset();
    else
    tabdata   = table();
    end
else
    data      = matrix{1};
    names     = matrix(2:end);
    if verLessThan('matlab', '8.3')
        tabdata   = dataset(matrix);  
    else
        tabdata   = array2table(data);
        tabdata.Properties.VariableNames = names;
    end
end