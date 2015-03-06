function [newdata, mask] =sliding_window_simple(data,window_size,thresh)
%sliding window algorithm


%get the row-by-row mean and SD
raw_mean = nanmean(data,2);
raw_sd = nanstd(data,0,2);

%quickly z-score each row by its mean and SD
temp_data = cellfun(@(x) (x-nanmean(x))/nanstd(x),num2cell(data,2),'Uniform',false);
temp_data = cell2mat(temp_data);



nx = size(temp_data,2);
nwind = window_size;

%this gives all the windows (window_size x width_of_data)
idx = bsxfun(@plus, (1:nwind)', 1+(0:(fix(nx/nwind)-1))*nwind)-1;

mask = zeros(size(temp_data));

%for each window
for i=1:size(temp_data,1)
    
    for k=1:size(idx,2)
        
        win=idx(:,k);
        
        if ~all(isnan(temp_data(i,win)))
            
            if abs(max(temp_data(i,win))-min(temp_data(i,win)))>thresh
                mask(i,win)=1;
            end
            
        end
        
    end
end

%remove the artifacts from the original data
newdata=data;
newdata(logical(mask)) = NaN;


end