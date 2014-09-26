function [clean_pupils,pupil_matrix] = pupil_remove_blinks(eyeStruct,varargin)


%make a copy of the original data
clean_pupils=eyeStruct;




for i=1:length(eyeStruct)
    
    %     fprintf('Trial %d:\n\n',i);
    
    blinks = eyeStruct(i).Blinks;
    
    
    if ~isempty(blinks)
        
        %         disp(length(blinks.sttime))
        
        for b=1:length(blinks.sttime)
            
            %get the beginning and end of blink
            blinkstart = max(blinks.sttime(b),1); %sometimes it's zero-- makes indexing crash
            blinkend = blinks.entime(b);
            
            
            %EyeLink will record a saccade starting just before, and ending
            %right after the blink. This is the real blink interval.
            saccades = eyeStruct(i).Saccades;
            
            
            sac_idx = [];
            
            if isfield(saccades,'sttime')
                %find the saccade that starts just before the blink
                sac_idx = find(saccades.sttime<=blinkstart,1,'last');
            end
            
            
            %if you can find it, use that for the timing
            if ~isempty(sac_idx)
                sac_start = saccades.sttime(sac_idx);
                sac_end = saccades.entime(sac_idx);
            else
                
                % blinks at the very beginning are a problem- saccades
                % aren't recorded. So, remove data from beginning until
                % first saccade. 
%                 if blinkstart<=100
%                     sac_start = 1;
%                     sac_end = saccades.sttime(1);
%                 else %if not at the beginning and somehow no saccade..
                    sac_start=max(blinkstart-25,1); %trim 25ms from beginning and end of blink
                    sac_end=min(blinkend+25,length(clean_pupils(i).pa));
%                 end
            end
            
            %remove that interval from the data
            clean_pupils(i).pa(sac_start:sac_end) = NaN;
            
          
        end
        
        %     else
        %         fprintf('No blinks found for trial %d',i);
        
    end
    
    
end

%sometimes it's useful to have everything in a flat matrix.

pupil_matrix={clean_pupils.pa}';
maxsamples=max(cellfun(@length,pupil_matrix));
pupil_matrix = cell2mat(cellfun(@(x) horzcat(x,nan(1,max(maxsamples-length(x),0))),pupil_matrix,'Uniform',false));


end


