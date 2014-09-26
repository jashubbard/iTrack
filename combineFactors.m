function [combined_factor, idx, varname] = combineFactors(allData,factors)



if ~isa(factors,'cell')
    factors=cellstr(factors);
end

nvars=length(factors);

allcols=cell(1,nvars);

%convert everything to strings
for i=1:nvars
    
    col=allData.(factors{i});
    
    if isa(col,'double') || isa(col,'logical')
        col=cellstr(num2str(col));
    end
    
    allcols{1,i}=col;
    
end


allcols=horzcat(allcols{:});

combined_factor=cell(size(allcols,1),1);

for i=1:size(allcols,1)
    
    temp=sprintf('%s_',allcols{i,:});
    temp=temp(1:end-1);
    combined_factor{i}=temp;
    
end

idx=grp2idx(combined_factor);

varname=sprintf('%s_',factors{:});
varname=varname(1:end-1);



end