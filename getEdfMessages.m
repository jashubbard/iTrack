function [onsets,itemOnset]=getEdfMessages(eyeStruct,varargin)


events=[eyeStruct.Events];
header=[eyeStruct.Header];



onsets = struct;
onsets.message=[events.message]'; %message like "STIMONSET"
onsets.time=[events.sttime]'-header.starttime;  %time that it occurred
%         ttrial=repmat(i,length(message),1); %add trial number

% onsets=dataset(message,time); %put it all in a dataset


emptyCells = cellfun(@isempty,onsets.message); %remove the empty rows
onsets.message(emptyCells,:) = [];
onsets.time(emptyCells,:) = [];

%         allonsets{i}=onsets;


if nargin>1
    msg=varargin{1};
    
    itemOnset=double(onsets.time(strcmp(onsets.message,msg)));
else
    itemOnset=NaN;
end




end