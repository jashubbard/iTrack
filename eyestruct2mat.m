function padata = eyestruct2mat(eyeStruct,varargin)


if ~isempty(varargin)
    field = varargin{1};
else
    field='pa';
end

%turn all the data into a matrix
temp={eyeStruct.(field)}';
maxsamples=max(cellfun(@length,temp));
padata=nan(length(eyeStruct),maxsamples);

for i=1:length(eyeStruct)
    padata(i,1:length(eyeStruct(i).(field))) = eyeStruct(i).(field);
end


end