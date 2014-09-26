function newCell = split_by_idx(mat,split_idx,varargin)
%splits matrix, mat, into a cell array based on unique values in split_idx
%if split_idx is a single number, it refers to a column of mat. If it's a
%vector, then mat will be split based on the unique values in that vector.
%The vector must have the same # of rows as mat.
%%

if length(split_idx)==1
    splitvec = mat(:,split_idx);
else
    splitvec = split_idx;
end

%take our vector of trial numbers to split the matrix
%into cells (1 for each trial)
% sortedA = vertcat(test{:});  %# Sort the rows by the first column
[~,~,uniqueIndex] = unique(splitvec);  %# Find indices of unique values
%#   in the first column
newCell = mat2cell(mat,accumarray(uniqueIndex(:),1),size(mat,2));  %#   into a cell array
% cellA{:}  %# Display the contents of the cells




%%
end