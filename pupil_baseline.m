function [eyeStruct,padata]= pupil_baseline(eyeStruct,baseline_times,varargin)

% padata=nan(length(eyeStruct),length(eyeStruct(1).pa_cleaned));



for i=1:length(eyeStruct)

pupildata=eyeStruct(i).pa;


baseline=nanmean(pupildata(baseline_times));
allsd = nanstd(pupildata(baseline_times(end):end));


%if entire baseline period is empty, OR the baseline mean is greater than
%20 SDs away from the rest of the data (from an artifact)
if isnan(baseline) || abs(baseline - min(length(eyeStruct(i).pa),nanmean(pupildata((baseline_times(end)+1):end)))) >= (2.0 * allsd) 
   
   pupildata(:) = NaN; %throw out everything
   
   %this attempts to fix things, but doesn't work under certain situations
%    firstsample=find(pupildata,1,'first');
%    baseline_window=baseline_times(end)-baseline_times(1);
%    baseline=nanmean(pupildata(firstsample:firstsample+baseline_window));
    
end

method='percent';

if nargin>2
    method=lower(varargin{1});
end

switch method
    case {'percent'}
        pupildata=(pupildata-baseline)./baseline;
    case {'subtract'}
        pupildata=(pupildata-baseline);
    case {'none'}
        pupildata=pupildata-0;
end


eyeStruct(i).pa_baselined=pupildata;
% padata(i,:)=pupildata;

end


padata=NaN;





end