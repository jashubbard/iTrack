function combined=combineDelimStrings(str1,str2,delim,varargin)
%interleaves 2 sets of delimited strings. Useful for making labels for
%legends, etc. 
%combined = combineDelimStrings('ID-Condition','123-2','-')
% returns 'ID-123-Condition-2';
%optional argument is a different delimiter, like "=":
%combined=combineDelimStrings('ID-'Condition','123-2','-','=')
%returns 'ID=123-Condition=2'
if nargin>3
    new_delim=varargin{1};
else
    new_delim=delim;
end


if isa(str1,'char')
    str1=cellstr(str1);
elseif ~isa(str1,'char') && ~isa(str1,'cell')
    str1=num2str(str1);
end


if isa(str2,'char')
    str2=cellstr(str2);
elseif ~isa(str2,'char') && ~isa(str2,'cell')
    str2=num2str(str2);
end


if any(size(str1)<size(str2))
    dimstochange=find(size(str1)<size(str2));
    
    s2size=size(str2);
    s2size(~dimstochange)=1;
    
    str1=repmat(str1,s2size);
    
elseif any(size(str2)<size(str1))
    dimstochange=find(size(str2)<size(str1));
    
    s1size=size(str1);
    s1size(~dimstochange)=1;
    
    str2=repmat(str2,s1size);
end



set1=cellfun(@(x) strsplit(x,delim),str1,'UniformOutput',false);
ncols1=max(cell2mat(cellfun(@length,set1,'Uniform',false)));


set2=cellfun(@(x) strsplit(x,delim),str2,'UniformOutput',false);
ncols2=max(cell2mat(cellfun(@length,set2,'Uniform',false)));


s1=cell(size(set1,1),ncols1);
s2=cell(size(set2,1),ncols2);

for i=1:size(set1,1)
    s1(i,:)=set1{i};
    s2(i,:)=set2{i};
end


temp_names = interleaveCols(s1,s2);

delim_row=cell(1,size(temp_names,2));

delim_row(1:2:end)={new_delim};
delim_row(2:2:end)={delim};

temp=interleaveCols(temp_names,repmat(delim_row,size(temp_names,1),1));
temp=temp(:,1:(end-1));

combined=cell(size(temp,1),1);


for i=1:size(temp,1)
   combined{i}=strjoin(temp(i,:),''); 
    
end

end