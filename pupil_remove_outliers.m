function clean=pupil_remove_outliers(padata,cutoff)

if isa(padata,'struct')
    
    dostruct = 1;
    eyeStruct = padata; %make a backup
    
    padata = eyestruct2mat(padata);
    
else
   dostruct = 0;
end

clean = nan(size(padata));


for row = 1:size(padata,1)
    data = padata(row,:);
    meandata=nanmean(data(:));
    sddata=nanstd(data(:));
    
    normdata=(data-meandata)/sddata;
    
    normdata(abs(normdata)>cutoff)=NaN;
    
    temp=(normdata.*sddata)+meandata;
    
    clean(row,:)=temp;
end

%if it's a structure, put everything back where it's supposed to go!
if dostruct==1
    
    for i=1:length(eyeStruct)
        eyeStruct(i).pa = clean(i,1:length(eyeStruct(i).pa));
    end
  
    clean = eyeStruct; %return a struct if a struct was given
end

end