 function dataDS=expand_table(tabledata,repeats)
    newvars = {'numfix','fix_idx'};
 
            %grab behavioral data
            tempbeh = table2cell(tabledata);
            
            %figure out how many fixations per trial
            numrepeats = repeats;
            numrepeats = repmat(numrepeats,1,size(tempbeh,2)+1);
            
            %add the column for "numrepeats" that gets replicated
            %for each trial (like all other behavioral data)
            tempbeh = horzcat(tempbeh,numrepeats(:,1));
            
            %create 1:num_fix for each trial
            fix_idx = cellfun(@(x) [1:x]',numrepeats(:,1),'Uniform',false);
            
            
            %replicate behavioral data for each fixation
            tempbeh = cellfun(@(x,y) repmat(x,y,1),tempbeh,numrepeats,'Uniform',false);
            
            %add the 1:num_fix column at the end
            tempbeh = horzcat(tempbeh,fix_idx);
            
            %turn into a table
            dataDS = cell2table(tempbeh,'VariableNames',horzcat(dataDS.Properties.VariableNames,newvars));
            
            
        end