function [clean,mask,cleanMat] = pupil_filter_blinks(padata,varargin)
%filters blinks based on algorithm focusing on differences between
%timepoints (i.e., not based on what EyeLink records). 

%INPUTS
%padata is a matrix where each row is a trial, or a structure from EyeLink with 'pa' as a field 
%optional argument (1/0) indicates whether to plot the raw and cleaned
%data. This is ignored if number of trials exceeds 25 (too many plots!)

%mask is a binary matrix/vector where 1 indicates that an artifact was found at that time


if nargin>1  %cutoff for number of plots
    doplot = varargin{1};
else
    doplot = 0;
end

if isa(padata,'struct')
    
    dostruct = 1;
    eyeStruct = padata; %make a backup
    
%     %turn all the data into a matrix
%     temp={eyeStruct.pa}'; 
%     maxsamples=max(cellfun(@length,temp));
%     padata=nan(length(eyeStruct),maxsamples);
%     
%     for i=1:length(eyeStruct)
%         padata(i,1:length(eyeStruct(i).pa)) = eyeStruct(i).pa;
%     end

    padata = eyestruct2mat(eyeStruct);

else
   dostruct = 0;
end

clean = nan(size(padata));
mask = zeros(size(clean));

numplots = 0;

%can be used for single trial or group of trials
for row = 1:size(padata,1)

enddata = max(find(~isnan(padata(row,:))))-1;
pa = padata(row,1:enddata);


% lb = 100;
%lower bound for determining artifacts 
%2 standard deviations below the mean
meanpa = nanmean(pa);
sdpa = nanstd(pa);

lb = meanpa - (2.5*sdpa);

clean(row,1:length(pa))= pa;
clean(row,pa<lb) = NaN;


dy1 = diff(pa);

dy1 = medfilt1(dy1,5);

[B,BN,BI] = RunLength(dy1); %downloaded from matlab file exchange. requires compiled MEX file
% [B,BN,BI] = RunLength_M(dy1); %this one doesn't require mex, a bit slower

blinks_ix = find(BN > (mean(BN)+3.0*std(BN)));

for j = 1:length(blinks_ix)
      
    blink_start = BI(blinks_ix(j));
    blink_end = BI(min(blinks_ix(j)+1,length(BI)));
    
    ix1 = blink_start:blink_end;
    
    %find decreasing values to the left of the blink
    i = find(dy1(1:min(ix1(2:end))) > 0,1,'last');
    ix2 = i:min(ix1(2:end));
    
    %remove the bad data!
    clean(row,ix2) = NaN;
    mask(row,ix2) =1;
    
    %find increasing values to the right of the blink
    i = find(dy1(max(ix1):(enddata-1)) < 0,1,'first');
    ix3 = max(ix1)+(0:i);
    
    %remove the bad data!
    clean(row,ix3) = NaN;
    mask(row,ix3) = NaN;

     
end

if doplot==1 && numplots <= 35
    figure;
    lastsample = find(~isnan(clean(row,:)),1,'last');
    x = 1:lastsample;
    plot(x,padata(row,1:lastsample),'--r',x,clean(row,1:lastsample),'g','LineWidth',2);
    legend({'raw','filtered'});
    numplots = numplots+1;

end

end

%if it's a structure, put everything back where it's supposed to go!
if dostruct==1
    for i=1:length(eyeStruct)
        eyeStruct(i).pa = clean(i,1:length(eyeStruct(i).pa));
    end
    
    cleanMat = clean; %save the matrix, it might be useful
    clean = eyeStruct; %but return a matrix as the first output
else
    cleanMat = []; %redundant if we don't have a structure
end

end