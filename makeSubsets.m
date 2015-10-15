function [subsets, levels, idx,varname]=makeSubsets(ds,factors,varargin)
%given an input table or dataset, splits into subsets based on factors,
%specified as a cell array of variable names
%optionally can specify only variables you want to keep in each subset
%and whether each subset should be another table, or just the data (i.e.,
%vectors or cell arrays, depending on the variable)
%subsets is a cell array of each subset. 
%levels is a table showing the unique groups (that match to subsets)
%idx is a vector of numbers as long as the original data specifying which
%row of "subsets" that row falls into 

p = inputParser;
p.addParameter('table',1,@(x) islogical(x) || ismember(x,[0,1]));
p.addParameter('keep',{},@(x) ischar(x) || iscell(x));
parse(p,varargin{:});

%only work with tables
if isa(ds,'dataset')
    ds = dataset2table(ds);
end

%which variables to keep
if isempty(p.Results.keep)
    varstokeep = ds.Properties.VariableNames;
else
    varstokeep = p.Results.keep;    
end


[levels,~,idx] = unique(ds(:,factors),'rows');

if p.Results.table
    subsets = arrayfun(@(x) ds(idx == x, varstokeep), unique(idx), 'Uniform', false);
else
    subsets = arrayfun(@(x) ds{idx == x, varstokeep}, unique(idx), 'Uniform', false);
end


levels.label = cell(height(levels),1);

[levels.label, varname,~]=makeLabels(levels,factors);



end