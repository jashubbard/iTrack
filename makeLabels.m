function [newcell, varnames,idx]=makeLabels(ds,vars,varargin)
%takes multiple cell arrays/vectors and combines them all into a single 
%cell array. Useful for combining multiple columns of a dataset into a
%single column.



nvars=length(vars);

allcols=cell(1,nvars);

if nargin==3
    varnames=varargin{1};
else
    varnames=vars;
end




if isa(varnames,'char')
    varnames=cellstr(varnames);
end

%convert everything to strings
for i=1:nvars
    
   col=ds.(vars{i});
   
   if isa(col,'double') || isa(col,'logical')
       col=cellstr(num2str(col));
   end
   
    col=addTextToVector(varnames{i},col,'=');
   
   allcols{i}=col;
    
end

allcols=horzcat(allcols{:});

newcell=cell(size(allcols,1),1);

for i=1:size(allcols,1)
    
    newcell{i}=combineWithDash(allcols(i,:)); %this slows everything down!
    
end


varnames=combineWithDash(varnames);

idx=grp2idx(newcell);

end

function newText=combineWithDash(strings)

newText=sprintf('%s-',strings{:});
newText=newText(1:end-1); %takes care of the Dash at the end


end