function [eyeStruct, artifact_mat] = pupil_sliding_window(eyeStruct,window_size,thresh)
%%
%given a matrix of pupil data, clean up artifacts by looking at
%peak-to-peak difference in a sliding window

%first, convert everything to a matrix to make it faster
data={eyeStruct.pa}';
maxsamples=max(cellfun(@length,data));
data = cell2mat(cellfun(@(x) horzcat(x,nan(1,max(maxsamples-length(x),0))),data,'Uniform',false));


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
% slidingWindows = temp_data(idx);

artifact_mat = ones(size(temp_data));

%for each window
for i=1:size(temp_data,1)

    for k=1:size(idx,2)
        
        win=idx(:,k);
        
        if ~all(isnan(temp_data(i,win)))
        
        if abs(max(temp_data(i,win))-min(temp_data(i,win)))>thresh
%             temp_data(i,win)=nanmean(temp_data(i,win));
%             temp_data(i,win)=NaN;
            artifact_mat(i,win)=NaN;
        end
        
        end
        %   cleandata=(temp_data.*nanstd(data))+nanmean(data);
        
    end
end


% cleandata = bsxfun(@times,temp_data,raw_sd);
% cleandata = bsxfun(@plus,cleandata,raw_mean); %bring it back into original units

%now remove all the artifacts from the data in eyeStruct
for i=1:length(eyeStruct)
   
    eyeStruct(i).pa = eyeStruct(i).pa.*(artifact_mat(i,1:length(eyeStruct(i).pa))); %multiply by the appropriate row in artifact_mat: it's all 1's, except NaN in positions where there are artifacts
    
end




%%
% close all
% randrows = randsample([150:size(cleandata,1)],50,0);
% 
% plot_all_rows(data(randrows,:))
% title('raw')
% 
% plot_all_rows(cleandata(randrows,:))
% title('cleaned')



end