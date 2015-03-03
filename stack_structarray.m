function [dataTable, dataStruct] = stack_structarray(s,varargin)
%take a structure array and "stack" each array so you have a scalar
%structure, then convert to a flat table. 
% d= struct;
% d(1).x = 5
% 


if nargin > 1
    fields = varargin{1};
else
    fields = fieldnames(s);
end

dataStruct = struct; 

 for f=1:length(fields)
       
        dataStruct.(fields{f}) = vertcat(s.(fields{f}));   
 end


dataTable = struct2table(dataStruct);




end


