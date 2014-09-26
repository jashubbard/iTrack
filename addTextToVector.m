function newCell = addTextToVector(left,right,varargin)
%combines text with some cell or vector. The order of inputs determines
%what goes on the left or right. Useful for generating labels for figure
%legends
% newcell = addTextToVector('x=',1:50) % returns a cell array where each
% item is "x=1","x=2" and so on. 
% newcell = addTextToVector('x=',1:50) %returns a cell array where each
% item is "1=x", "2=x", and so on. 


if nargin==3
    delim=varargin{1};
else
    delim='-';
end

if isnumeric(left)  
    left = cellfun(@num2str,num2cell(left),'UniformOutput',false);
elseif ischar(left)
    left = cellstr(left);
end

if isnumeric(right)
    
    right = cellfun(@num2str,num2cell(right),'UniformOutput',false);
elseif ischar(right)
    right = cellstr(right);
end


if any(size(left)<size(right))
    dimstochange=find(size(left)<size(right));
    
    rsize=size(right);
    rsize(~dimstochange)=1;
    
    
    
    
    left=repmat(left,rsize);
    
elseif any(size(right)<size(left))
    dimstochange=find(size(right)<size(left));
    
    lsize=size(left);
    lsize(~dimstochange)=1;
    
    right=repmat(right,lsize);
end

if size(left,1)==1
    left=left';
end

if size(right,1)==1
    right=right';
end

temp=horzcat(left,right);

newCell=cell(size(temp,1),1);

for i=1:size(temp,1)
   
    newCell{i}=combineWithDelim(temp(i,:),delim);
    
    
end




end


function newText=combineWithDelim(strings,delim)


newText=sprintf(['%s' delim],strings{:});
newText=newText(1:end-1); %takes care of the Dash at the end



end