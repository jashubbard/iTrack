function [allsubsets,idx,levels]=makeSubsets_OLD(allData,factors,varargin)

if nargin>2
   
    names=varargin{1};
    
else
    names=factors;
    
end


[combined_factor,idx,varname] = combineFactors(allData,factors);

levels=unique(idx);
levels(isnan(levels))=[];

allsubsets=struct;

allData.(varname)=combined_factor;

for i=1:length(levels)
    allsubsets(i).subset=allData(idx==levels(i),:);
    val=unique(combined_factor(idx==levels(i)));
    allsubsets(i).group=val{1}; 
    allsubsets(i).varname=varname;
    allsubsets(i).grouping_vars=factors;
    allsubsets(i).title = combineDelimStrings(varname,val{1},'_','=');
end



end