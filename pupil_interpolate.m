function clean = pupil_interpolate(padata,varargin)
%filters blinks based on algorithm focusing on differences between
%timepoints (i.e., not based on what EyeLink records). 

%INPUTS
%padata is a matrix where each row is a trial, or a structure from EyeLink with 'pa' as a field 
%optional argument (1/0) indicates whether to plot the raw and cleaned
%data. This is ignored if number of trials exceeds 25 (too many plots!)

%mask is a binary matrix/vector where 1 indicates that an artifact was found at that time
warning('off','MATLAB:chckxy:IgnoreNaN');

%can specify whether the smoothed data replaces everything, or just fills
%in the gaps (default)
if nargin>1
    replace = varargin{1};
else
    replace = 0;
end



if isa(padata,'struct')
    
    dostruct = 1;
    eyeStruct = padata; %make a backup
    
    padata = eyestruct2mat(padata);
%     %turn all the data into a matrix
%     temp={eyeStruct.pa}'; 
%     maxsamples=max(cellfun(@length,temp));
%     padata=nan(length(eyeStruct),maxsamples);
%     
%     for i=1:length(eyeStruct)
%         padata(i,1:length(eyeStruct(i).pa)) = eyeStruct(i).pa;
%     end
    
else
   dostruct = 0;
end

clean = nan(size(padata));


%can be used for single trial or group of trials
for row = 1:size(padata,1)
    
    startdata = find(~isnan(padata(row,:)),1,'first');
    enddata = find(~isnan(padata(row,:)),1,'last');
    pa = padata(row,1:enddata);
    
    
    if mean(~isnan(padata(row,startdata:enddata))) >=.30 && length(startdata:enddata)>5
        temp = pa(startdata:enddata);
        sm = splinefit(1:length(temp),temp,25);
        sm = ppval(sm,1:length(temp));
        
        if replace
            clean(row,startdata:enddata)=sm;
        else
            
            temp(isnan(temp)) = sm(isnan(temp));
            clean(row,startdata:enddata) = temp;
        end
        
    else
        clean(row,1:length(pa)) = pa;
    end
    

end

%if it's a structure, put everything back where it's supposed to go!
if dostruct==1
    for i=1:length(eyeStruct)
        eyeStruct(i).pa = clean(i,1:length(eyeStruct(i).pa));
    end
    
%     cleanMat = clean; %save the matrix, it might be useful
    clean = eyeStruct; %but return a matrix as the first output
% else
%     cleanMat = []; %redundant if we don't have a structure
end

warning('on','all');
end