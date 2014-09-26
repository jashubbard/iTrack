function avdata=moving_average(y,window_width,varargin)

if ~isempty(varargin)
    type = varargin{1};
else
    type = 'same';
end
%%
nan_idx = isnan(y);
y(nan_idx) = 0;

mask=ones(1,window_width)/window_width;

avdata=conv2(y,mask,type);

%the ends don't work correctly
avdata(1:window_width)=nanmean(y(1:window_width));
avdata(end-window_width:end)=nanmean(y(end-window_width:end));
%%

if strcmp(type,'same')
    avdata(nan_idx) = NaN;
end


end