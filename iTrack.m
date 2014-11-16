classdef iTrack
    
    properties
%         fields={'subjects','beh','eyedata'};
        subs=cell(1,1);
        data=cell(1,1);
        raw=cell(1,1);
        subject_var='';
        screen = struct;
        rois = struct;
        log = struct;
        
        
    end
    
    methods
        
        function obj=iTrack(varargin)
            p = inputParser;
            p.addParameter('edfs',{},@(x) iscell(x) || ischar(x));
            p.addParameter('samples',0,@(x) islogical(x) || ismember(x,[0,1]))
            p.addParameter('subject_var','ID',@ischar)
            p.addParameter('resolution',[1024,768],@(x) length(x)==2 && isnumeric(x));
            p.addParameter('overlay',nan(1,3),@(x) isnumeric(x) && size(x,1)==3);
            p.addParameter('keepraw',false);
            p.addParameter('keeptrials',[]);
            p.addParameter('droptrials',[]);
            
            parse(p,varargin{:});
            
            edfs = p.Results.edfs;
            
            obj.subject_var = p.Results.subject_var;
            
            obj.screen.dims=p.Results.resolution;
            
            obj.rois.single = struct;
            obj.rois.combined = struct;
            
            
                %1 argument, cell of file names including full path
            if ~isempty(edfs)
                path=pwd;
                             
               
                
            else
                %default, no arguments, prompt for edf files
                [edfs, path] = uigetfile({'*.edf'},'select edf files','MultiSelect','on');
            end

            cd(path);
            
            %if only 1 filename given
            if ~isa(edfs,'cell')
                edfs={edfs};
            end
            
            %main place where data is kept
            obj.data=cell(length(edfs),1);
            
            %loop through all edfs
            for s=1:length(edfs)
                
                fname=edfs{s};
                
                
                %import edf file
                
                if p.Results.samples==1
                    x=edfImport(fname,[1 1 1],'pa gx gy');
                else
                    x=edfImport(fname,[1 1 0]);
                end
                
                x=edfExtractInterestingEvents(x);
                
                
                eye=x(1).Header.rec.eye;
                sample_rate=x(1).Header.rec.sample_rate;
                
                %if we want to throw out some of the data
                %for now, the same for each subject
                
                %keeping only certain trials
                if ~isempty(p.Results.keeptrials)
                    x = x(p.Results.keeptrials);
                end
                
                %dropping certain trials
                if ~isempty(p.Results.droptrials)
                    x(p.Results.droptrials)=[];
                end
                
                %add these fields even if they're not necessarily used
                [x.sample_rate] = deal(sample_rate);
                [x.pa] = deal([]);
                [x.gx] = deal([]);
                [x.gy] = deal([]);
                [x.fixation_times] = deal([]);
                [x.fixation_time_mask] = deal([]);
                [x.beh] = deal(struct);
%                 [x.fixation_hits]  = deal([]);
                
                %trials with no recorded fixations are empty [], but we want
                %them to be structures with NaN's for each field - it makes
                %life easier

                %create a blank structure for trials with no fixations
                fnames = fieldnames([x.Fixations]); %fieldnames from fixation structure
                
                empty_struct = struct;
                
                for f=1:length(fnames)
                    empty_struct.(fnames{f}) = NaN;
                end
                
               
                %fill those trials with our empty struct
                [x(cellfun(@isempty,{x.Fixations})).Fixations] = deal(empty_struct);
               
                %this transposes all the fixation data so it's 1 row per
                %fixation, instead of column-wise
                temp = arrayfun(@(y) structfun(@transpose,y,'Uniform',false),[x.Fixations],'Uniform',false);
                [x.Fixations] = deal(temp{:});
                
                for i=1:length(x)
                    
                    x(i).events=getEdfMessages(x(i));%external file
                    
                    if p.Results.samples
                        numsamples = length(x(i).Samples.pa); %when importing samples, a few extra samples are trimmed off the end
                    else
                        numsamples = x(i).Header.duration; %this is accurate enough for fixation data
                    end
                    
                    x(i).index.raw=single(1:numsamples);
                    x(i).numsamples = numsamples;
                    
                    
                    %create a matrix of zeros for each trial that is 
                    % number fixations x  number of samples. Used for indexing fixation data
                    if ~isempty(x(i).Fixations)
                       x(i).fixation_times =single(zeros(length(x(i).Fixations.sttime),x(i).numsamples));
                    end
                    
                    if p.Results.samples==1
%                         temp=x(i).Samples;
%                         x(i).index.raw=double(0:length(temp.pa(eye,:))-1);
%                         x(i).numsamples=x(i).index.raw(end)+1;
                        x(i).pa=x(i).Samples.pa(eye,:);
                        x(i).gx=x(i).Samples.gx(eye,:);
                        x(i).gy=x(i).Samples.gy(eye,:);
                        
                    end
                end
                
                if p.Results.samples==1
                    x = rmfield(x,'Samples'); %removing redundant data
                end
                
                x= rmfield(x,'Events'); %also remove this to save memory
                x= rmfield(x,'Buttons'); %also remove this to save memory
                
                obj.raw{s,1}=x; %structure with all the information
                obj.subs{s,1}=str2double(regexprep(fname,'[!@#$%^&()?"*+=-_'',./~` A-Za-z]','')); %subject number from filename, minus any funny characters
                
                clear x;
            end
            
            %copy everything from the raw structure into our "data"
            %structure. That way we can always revert back to the raw data
            %without reloading edfs. 
            obj.data=obj.raw; 
            
            if ~p.Results.keepraw
                obj.raw = [];
            end
            
            %check to see if everyone was recorded at same sample rate
            checksamprate = unique(cellfun(@(x) unique(cell2mat(x)),subsref(obj,struct('type','.','subs','sample_rate')),'Uniform',true));
            if length(checksamprate)>1
                warning('WARNING: You have different sample rates for different subjects! Use "resample" if you are looking at pupil data!');
            end
            
            
            obj = index_fixations(obj); %represents the fixations with binary matricies, useful for other operations
            
            disp('done!');
        end
        
        
        function obj=reset(obj)
         %revert back to the raw data without reloading edfs. 
            obj.data=obj.raw;
            obj = index_fixations(obj);

        end

        function obj= add_behdata(obj,behFile,varargin)
            p = inputParser;
            p.addParameter('index',{},@iscell);
            p.addParameter('addevent',false);
            p.addParameter('append',false);
            parse(p,varargin{:});
            
            
            numsubs=size(obj.data(:,1),1);
            
            
            if ~isempty(p.Results.index)
               index_vars = p.Results.index; 
            else
               index_vars = {};
            end
            
            
            
            
            %can handle files (txt,dat, or csv) or variables (dataset or
            %table)
            
            if isa(behFile,'table')
                 behDS = behFile;
            elseif isa(behFile,'dataset')
                 behDS = dataset2table(behFile);
            
            elseif exist(behFile,'file')
                [~,~,ext] = fileparts(behFile);
                
                switch ext
                    case {'.txt','.dat'} %assumes tab-delimited with header row
                        behDS = readtable(behFile,'Delimiter','\t');
                    case {'.csv'} %assumes comma-delimited with header row
                        behDS = readtable(behFile,'Delimiter',',');
                end
           
            end
            
            
            
            if p.Results.append
                [~,existing_vars] = get_all_fields(obj);
                newvars = setdiff(behDS.Properties.VariableNames,existing_vars);
            
            end

            
            
            for s=1:numsubs
                
                eyeStruct=obj.data{s};
                eyeDS=table;
                
                behSubset=behDS(behDS.(obj.subject_var)==obj.subs{s},:);
                
                %subject variable (ID/SUBID) is specified earlier.
                %add a column for that subject
                eyeDS.(obj.subject_var)=repmat(obj.subs{s},length(eyeStruct),1);
                eyeDS.eye_idx = [1:height(eyeDS)]'; %add the eye_idx variable in case there is eyetracking data that doesn't match with behavioral
                
                
                if ~isempty(index_vars)
                    %if index vars specified
                    
                    %copy from eyeStruct into eyeDS
                    for i=1:length(index_vars)
                        
                        cur_var=index_vars{i};
                        var=[eyeStruct.(cur_var)]';
                        
                        eyeDS.(cur_var)=var;
                        
                    end
                    

                    
                    joinedDS=outerjoin(eyeDS,behSubset,'Keys',horzcat(obj.subject_var,index_vars),'MergeKeys',true);
                else
                    
                    behSubset.eye_idx = [1:height(behSubset)]';
                    
                    if height(eyeDS) ~= height(behSubset)
                        error('Index variables not specified, and number of rows do not match the eye data!');
                    end
                    
                    joinedDS=outerjoin(eyeDS,behSubset,'Keys',horzcat(obj.subject_var,{'eye_idx'}),'MergeKeys',true);
                    
                end
                
                
                
                joinedDS = sortrows(joinedDS,'eye_idx');                
                joinedDS = table2struct(joinedDS);
                
                for i = 1:length(eyeStruct)
                    
                    if p.Results.append
                        
                        for v = 1:length(newvars)
                            
                            %for adding events, instead of behavioral data
                            if p.Results.addevent
                                
                                if ~isnan(joinedDS(i).(newvars{v})) %if missing data, don't add the event!
                                eyeStruct(i).events.message{end+1}=newvars{v};
                                eyeStruct(i).events.time(end+1)=joinedDS(i).(newvars{v});
                                end
                                
                            else    
                            
                            eyeStruct(i).beh.(newvars{v}) = joinedDS(i).(newvars{v});
                            
                            end
                        end
                    else
                        eyeStruct(i).beh = joinedDS(i);
                    end
                end     
                    
                obj.data{s}=eyeStruct;
                
                
            end
            
        end
        
        
        
        
        function obj=resample(obj,desired_samplerate)
            
            numsubs=size(obj.data(:,1),1);
            
            
            for s=1:numsubs
                eyeStruct=obj.data{s};
                
                
                sample_rate=eyeStruct(1).sample_rate;
                
                if sample_rate ~= desired_samplerate
                    fprintf('Resampling subject %d from %d to %d Hz\n',s,sample_rate,desired_samplerate);
                    
                    
                    if sample_rate< desired_samplerate
                        interp_direction='up';
                        interp_factor= desired_samplerate/sample_rate;
                    elseif sample_rate>desired_samplerate
                        interp_direction='down';
                        interp_factor=sample_rate/desired_samplerate;
                        
                    end
                    
                    
                    
                    for i=1:length(eyeStruct)
                        eyeStruct(i).rawsample.pa=eyeStruct(i).pa;
                        eyeStruct(i).rawsample.gx=eyeStruct(i).gx;
                        eyeStruct(i).rawsample.gy=eyeStruct(i).gy;
                        
                        switch interp_direction
                            case {'up'}
                                eyeStruct(i).pa=interp(eyeStruct(i).pa,interp_factor);
                                eyeStruct(i).gx=interp(eyeStruct(i).gx,interp_factor);
                                eyeStruct(i).gy=interp(eyeStruct(i).gy,interp_factor);
%                                 eyeStruct(i).onsets.time = fix(eyeStruct(i).onsets.time.*interp_factor);
                            case {'down'}
                                eyeStruct(i).pa=downsample(eyeStruct(i).pa,interp_factor);
                                eyeStruct(i).gx=downsample(eyeStruct(i).gx,interp_factor);
                                eyeStruct(i).gy=downsample(eyeStruct(i).gy,interp_factor);
%                                 eyeStruct(i).onsets.time = fix(eyeStruct(i).onsets.time.*interp_factor);
                        end
                        
                        
                        eyeStruct(i).sample_rate=desired_samplerate;
                        eyeStruct(i).numsamples=length(eyeStruct(i).pa);
                        
                        
                    end
                    
                    
                    fprintf('Subject %d complete\n',s);
                
                else
                    
                    for i=1:length(eyeStruct)
                        
                    eyeStruct(i).rawsample.pa=eyeStruct(i).pa;
                    eyeStruct(i).rawsample.gx=eyeStruct(i).gx;
                    eyeStruct(i).rawsample.gy=eyeStruct(i).gy;
                    
                    end
                
                end
                
                obj.data{s}=eyeStruct;
                
            end
            
            
            
            
        end
        
        
        
     
        
        function obj=baseline(obj,baseline_times,varargin)
             p = inputParser;
            p.addParameter('method','percent',@ischar);
            p.addParameter('func','mean',@(x) ismember(x,{'mean','mad'}));
             p.addParameter('raw',true,@(x) islogical(x) || ismember(x,[0,1]));
            p.addParameter('event','',@ischar);
            parse(p,varargin{:});
            
            

            
            
            
            numsubs=size(obj.data(:,1));
            

            
            for s=1:numsubs
                
                
                
                
                eyeStruct=obj.data{s};
                
                if eyeStruct(1).sample_rate ~=1000
                    
                    sr_fraction = eyeStruct(1).sample_rate/1000;
                    
                    baseline_interval = [baseline_times(1)*sr_fraction:baseline_times(end)*sr_fraction]; 
                else
                    baseline_interval = baseline_times;
                    
                end
                
                
                for i=1:length(eyeStruct)
                    
                    if ~p.Results.raw
                        pupildata = eyeStruct(i).aligned.(p.Results.event);
                    else
                        
                        pupildata=eyeStruct(i).pa;
                    end
                    
                    
                    switch(p.Results.func)
                        case {'mean'}
                            baseline=nanmean(pupildata(baseline_interval));
                        case {'mad'}
                            baseline=mad(pupildata(baseline_interval));
                    end

%                     allsd = nanstd(pupildata(baseline_times(end):end));
                    
                    numsamples = length(baseline_interval);
                    
                    %if entire baseline period is empty, OR the baseline mean is greater than
                    %20 SDs away from the rest of the data (from an artifact)
    %                     if isnan(baseline) || abs(baseline - nanmean(pupildata(min(length(eyeStruct(i).pa),baseline_times(end)+1):end))) >= (2.0 * allsd)
                    if isnan(baseline)
%                         pupildata(:) = NaN; %throw out everything
                        
                        %this attempts to fix things, but doesn't work under certain situations
                           firstsample=find(pupildata,1,'first');
                           baseline = nanmean(pupildata(firstsample:(firstsample+numsamples-1)));
                        
                    end
                    
                    
                    switch p.Results.method
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
                
                
%                 padata=NaN;
                
                
                
                
                obj.data{s}=eyeStruct;                
                
            end
            
        end
        
        function obj=index_events(obj,varnames,varargin)
            
            
            numsubs=size(obj.data(:,1));
            
            if ~isa(varnames,'cell')
                varnames={varnames};
            end
            
            
            for s=1:numsubs
                
                eyeStruct=obj.data{s};
                
                for t=1:length(eyeStruct)
                    
                    trial=eyeStruct(t);
                    
                    for v=1:length(varnames)
                        
                        var1=varnames{v};
                        
                        
                        vartime=double(trial.events.time(find(ismember(trial.events.message,var1))));
                        
                        if ~isempty(vartime)
                            start=double(1-vartime);
                            endtime=double(trial.numsamples-vartime);
                            
                            eyeStruct(t).index.(var1)=single(start:endtime);
                        else
                            eyeStruct(t).index.(var1)=single(nan(1,trial.numsamples));
                            fprintf('Event %s not found for subject %d, trial %d\n',var1,s,t);
                        end
                        
                    end
                    
                end
                
                obj.data{s}=eyeStruct;
                
            end
            
            
        end
        
        function obj=lineup(obj,lineup_var,before_lineup,after_lineup,varargin)
            %to make old scripts work. 
            obj = epoch(obj,lineup_var,before_lineup,after_lineup,varargin{:});
            
        end
        
        
        function obj=epoch(obj,varargin)
            p = inputParser;
            p.addParameter('type','baselined',@(x) ismember(x,{'fix','baselined','raw'}));
            p.addParameter('event','',@ischar);
            p.addParameter('rois',{'none'},@iscell);
            p.addParameter('interval',[],@isnumeric);
            
            parse(p,varargin{:});
            
            before_lineup = p.Results.interval(1);
            after_lineup = p.Results.interval(2);
            lineup_var = p.Results.event;
            
            
           
            
            
            numsubs=length(obj.subs);
            
            
            for s=1:numsubs
                
                eyeStruct=obj.data{s};
                
             for r = 1:length(p.Results.rois)   
                for t=1:length(eyeStruct)
                    
                    trial=eyeStruct(t);
                    
                    
                    if ~all(isnan(trial.index.(lineup_var)))
                        
                        %                         startwindow=find(window(1)==trial.index.(lineup_var));
                        %                         endwindow=find(window(end)==trial.index.(lineup_var));
                        %%
                        startwindow=find(trial.index.(lineup_var)==((before_lineup-1)*-1)); %time of event is 0, so 500 ms before would be coded as -499 (NOT -500), hence the -1
                        endwindow=find(trial.index.(lineup_var)==after_lineup);
                        
%%
                        
                        switch p.Results.type
                            case {'baselined'}
                                newdata=trial.pa_baselined;
                            case {'raw'}
                                newdata=trial.pa;
                            case {'fix'}
                                newdata=trial.fix;
                                newdata = newdata.(p.Results.rois{r});
                        end
                        
                        if isempty(startwindow)
                            difference=abs(double(trial.index.(lineup_var)(1))-((before_lineup-1)*-1));
                            newdata=horzcat(nan(1,difference),newdata);
                            startwindow=1;
                            endwindow=endwindow+difference;
                        end
                        
                        
                        if isempty(endwindow)
                            %                             difference=abs(window(end)-varidx(end));
                            difference=abs(double(trial.index.(lineup_var)(end))-after_lineup);
                            newdata=horzcat(newdata,nan(1,difference));
                            endwindow=length(newdata); %-1 here?
                        end
                        
                        
                        newdata=newdata(startwindow:endwindow);
                        
                        if length(newdata) ~= length((before_lineup-1)*-1:(after_lineup));
                            fprintf('Warning: subject %d, trial %d has length of %d\n',obj.subs{s},t,length(newdata));
                        end
                        
                    else %if event isn't found, just write NaNs
                        
                        window=((before_lineup-1)*-1):(after_lineup);
                        window_width=length(window);
                        newdata=nan(1,window_width);
                    end
                    
                    if ~strcmp(p.Results.type,'fix')
                        eyeStruct(t).windows.(lineup_var)=[before_lineup*-1, after_lineup];
                        eyeStruct(t).aligned.(lineup_var)=newdata;
                    else
                        eyeStruct(t).fix_windows.(lineup_var)=[before_lineup*-1, after_lineup];
                        eyeStruct(t).fix_aligned.(lineup_var).(p.Results.rois{r})=newdata;
                    
                    end
                    
                    
                end
                
            end
                obj.data{s}=eyeStruct;
                
                
                
            end
            
            
            
        end
        
        
        function allfigdata = linePlots(obj,varargin)
           %wrapper function for linePlots (which can be used
           %independently)
           %give it variables for separating plots (plotVars) and for
           %separating lines (lineVars), and the data that you want (dataVar)
            p = inputParser;
            p.addParameter('plotVars',{},@iscell);
            p.addParameter('lineVars',{},@iscell);
            p.addParameter('fix',true,@(x) islogical(x) || ismember(x,[0,1]));
            p.addParameter('dataVar','',@ischar);
            p.addParameter('event','',@ischar);
            p.addParameter('roi','',@ischar);
            p.addParameter('legend',true,@(x) islogical(x) || ismember(x,[0,1]))
            p.addParameter('LineWidth',2,@isnumeric)
            p.addParameter('ylim','auto',@(x) (length(x)==2 && isnumeric(x)) || strcmp(x,'auto'));
            p.addParameter('ttest',false,@(x) islogical(x) || ismember(x,[0,1]))
            p.addParameter('alpha',.05,@isnumeric)
            p.addParameter('filter','',@ischar);
            p.addParameter('labels',{},@iscell);
            p.addParameter('func',@nanmean);
            p.addParameter('downsample',[],@isnumeric)
            p.addParameter('value',1);
            
            parse(p,varargin{:});
            
            subVar = obj.subject_var;
            plotVars = p.Results.plotVars;
            lineVars = p.Results.lineVars;
            dataVar = p.Results.dataVar;
            
            if ~isempty(p.Results.lineVars)
                allfactors = horzcat(subVar,plotVars,lineVars);
            else
                allfactors = horzcat(subVar,plotVars);
            end
            
            if ~isempty(p.Results.filter)
                allfactors = horzcat(allfactors,p.Results.filter);
            end
            
            
            if p.Results.fix
                allData = build_dataset(obj,'beh',allfactors,'events',{dataVar},'fix',true,'rois',{p.Results.roi});
                dataVar = strcat(p.Results.roi,'_',p.Results.dataVar);
                
            else
                allData = build_dataset(obj,'beh',allfactors,'events',{dataVar});
            end
            
            
            if ~isempty(p.Results.filter)
                allData = allData(allData.(p.Results.filter)==p.Results.value,:);
            end
            
            
            if ~isempty(p.Results.downsample)
                allData.(dataVar) = downsample(allData.(dataVar)',p.Results.downsample)';
            end
            
            
            allfigdata = linePlots(allData,'plotVars',p.Results.plotVars,...
                'lineVars',p.Results.lineVars,...
                'dataVar',dataVar,...
                'subVar',subVar,...
                'fix',p.Results.fix,...
                'roi',p.Results.roi,...
                'legend',p.Results.legend,...
                'ylim',p.Results.ylim,...
                'LineWidth',p.Results.LineWidth,...
                'ttest',p.Results.ttest,...
                'alpha',p.Results.alpha,...
                'func',p.Results.func);
            
            
        end
        
        
        
        function [data,datamat]=get(obj,varargin)
           %generic function for pulling data from iTrack objects
           %options are:
           %'aligned' - epoched data, lined up to some event
           %'raw' - raw pupil data
           %'allfix' - data for every fixation
           %'messages' - all events that are sent to eyelink during the experiment
           %
           % returns 2 outputs- first is a N x 2 cell array for every
           % subject. First column is a vector of subject's ID, second
           % column is the requested data for that subject. 
           % in some cases it also returns a matrix of the requested data,
           % with the first column being the subject ID. 
           %
           % [dataCell, dataMat] = get(obj,'aligned','event','STIMONSET');
           % [dataCell, dataMat] = get(obj,'raw');
           
            p = inputParser;
            p.addRequired('item');
%             p.addRequired('eyevars');
            p.addParameter('event',[]);
            p.addParameter('behvars',{},@iscell);
            p.addParameter('search','');
            p.addParameter('roi','all');
            p.addParameter('fix',false);
            parse(p,varargin{:});
            
            
            item = p.Results.item;
            
            numsubs=size(obj.data(:,1),1);
            
            data=cell(numsubs,2);
            
            %if we're getting all fixation data, this only needs to be
            %called once. 
            if ~isempty(p.Results.roi)
                if strcmp(p.Results.roi,'xy')
                    temp_allfix = get_Fixations(obj,'fields',{'gavx','gavy'});
                    
                    
                elseif strcmp(p.Results.roi,'time')
                    temp_allfix = get_Fixations(obj,'fields',{'sttime','entime'});
                end
            end
            
            
            for s=1:numsubs
                
                switch item
                    case {'aligned'}
                        lineup_var=p.Results.event;
                        
                        temp=obj.data{s};
                        
                        if p.Results.fix==false
                            temp=[temp.aligned];
                            temp=vertcat(temp(:).(lineup_var));
                        else
                            temp = [temp.fix_aligned];
                            temp = {temp.(p.Results.event)};
                            temp = cellfun(@(x) x.(p.Results.roi),temp,'Uniform',false)';
                            temp = vertcat(temp{:});
                        end
                        
                        
                    case {'raw'}
                        padata=obj.data{s};
                        padata={padata.pa}';
                        maxsamples=max(cellfun(@length,padata));
                        
                        temp=nan(size(padata,1),maxsamples);
                        
                        
                        for i=1:size(padata,1)
                            
                            temp(i,1:length(padata{i}))=padata{i};
                            
                            
                        end
                    case {'events'}
                                                
                        all_messages = arrayfun(@(x) x.events.message',obj.data{s},'Uniform',false);
                        all_times = arrayfun(@(x) x.events.time',obj.data{s},'Uniform',false);
%                         temp = horzcat(all_messages{:})'; %for stacking
%                         all messages with no organization
                        temp = all_messages';
                        all_times = all_times';


                        if nargin>2
                            search_exp = p.Results.search;
                            
%                             all_matches = ~cellfun(@isempty,regexp(temp,search_exp));
                            
                            temp2 = cell(size(temp,1),2);
                            
                            for i=1:length(temp2)
                                
                                msg = temp{i};
                                msgtime = all_times{i};
                                
                                all_matches = ~cellfun(@isempty,regexp(msg,search_exp));
                                
                                temp2{i,1} = msg(all_matches);
                                temp2{i,2} = msgtime(all_matches);
                            end
                            
                            
                            temp = temp2;
                        end
                      case {'allfix'}
                       
                          %getting the xy data or timing data is a special
                          %case
                          
                          if strcmp(p.Results.roi,'xy') || strcmp(p.Results.roi,'time')
                              temp = temp_allfix{s};
                              
                          else
                              temp= obj.data{s};
                              temp = [temp.fixation_hits];
                              temp = reshape({temp.(p.Results.roi)},[],1);
                          end
                          
                end
                
                
                
                data{s,1}=repmat(obj.subs{s},size(temp,1),1);
                data{s,2}=temp;
                
            end
            
            
            switch item
                case {'raw'}
                    
                    ncols=max(cellfun(@(x) size(x,2),data(:,2)));
                    nrows=max(cellfun(@(x) size(x,1),data(:,2)));
                    
                    datamat=nan(nrows,ncols);
                    
                    for i=1:size(data,1)
                       difference= ncols-size(data{i,2},2);
                       
                       padding=nan(size(data{i,2},1),difference);
                       
                       data{i,2} = horzcat(data{i,2},padding);
                        
                    end
                    
                    datamat=cell2mat(vertcat(data(:,:)));
                    

                case {'messages'}
                    
                    all_subs = vertcat(data{:,1});
                    all_msg = cellfun(@cell2mat,vertcat(data{:,2}),'Uniform',false);
                    
                    datamat = nan(length(all_subs),2);
                    datamat(:,1) = all_subs;
                    
                    
                    for i =1:length(all_msg)
                        
                        msg = all_msg{i};
                        
                        if isa(msg,'cell') && ~isempty(msg)
                            msg = msg{1}; %only grab the first one;
                        end
                        
                        
                        if isa(msg,'char') && ~isempty(msg)
                         
                            num = str2double(msg(regexp(msg,'[\d.]')));
                            
                            if ~isempty(num)
                                datamat(i,2) = num;
                            end
                        elseif isa(msg,'double') && ~isempty(msg)
                            datamat(i,2) = msg;
                        end
                            
                    end
                case {'allfix'}
                    datamat = [];
                        
                otherwise
                    datamat=cell2mat(vertcat(data(:,:)));
            end
            
          
            
            
        end
        
        
        
        function dataDS=build_dataset(obj,varargin)
            
            p = inputParser;
            p.addParameter('fix',false,@(x) islogical(x) || ismember(x,[0,1]));
            p.addParameter('allfix',false,@(x) islogical(x) || ismember(x,[0,1]));
            p.addParameter('raw',false,@(x) islogical(x) || ismember(x,[0,1]));
            p.addParameter('beh',{},@iscell);
            p.addParameter('events',{'none'},@iscell);
            p.addParameter('rois',{},@iscell);
            p.addParameter('downsample',[],@isnumeric)
            parse(p,varargin{:});
            
            
            eyevars = p.Results.events;
            
            if ~isa(eyevars,'cell')
                eyevars={eyevars};
            end
            

            %             for s=1:numsubs
            
            allsubdata=[obj.data{:}]; %behavioral data
            behdata=struct2table(vertcat(allsubdata.beh));
            
            %                 behdata=behdata(:,behvars);
            
            %by default, grab all behavioral data. Otherwise, specify
            %variables
            if isempty(p.Results.beh)
                dataDS = behdata;
            else
                dataDS=behdata(:,p.Results.beh);
            end
            
            
            for v=1:length(eyevars)
                
                if p.Results.raw
                        %for grabbing raw pupil data
                        [~, temp] = get(obj,'raw');
                        fname = eyevars{v};
                    
                elseif p.Results.allfix
                        %for grabbing all fixation data (multiple fixations
                        %per trial
                        temptable = table;
                        
                        for r=1:length(p.Results.rois)
                            
                            if ~strcmp(p.Results.rois{r},'xy') && ~strcmp(p.Results.rois{r},'time')
                                roiname = strcat('roi_',p.Results.rois{r});
                                suffix = '_hit';
                            else
                                roiname = p.Results.rois{r};
                                suffix = '';
                            end
                            
                            temp = get(obj,'allfix','roi',roiname);
                            fname = strcat(strrep(p.Results.rois{r},'roi_',''),suffix);
                            
                            if strcmp(p.Results.rois{r},'xy')
                                temptable.x_coord = cellfun(@(x) x(:,1),vertcat(temp{:,2}),'Uniform',false);
                                temptable.y_coord = cellfun(@(x) x(:,2),vertcat(temp{:,2}),'Uniform',false);
                            elseif strcmp(p.Results.rois{r},'time')
                                temptable.fixstart = cellfun(@(x) x(:,1),vertcat(temp{:,2}),'Uniform',false);
                                temptable.fixend = cellfun(@(x) x(:,2),vertcat(temp{:,2}),'Uniform',false);
                            else
                                temptable.(fname) = vertcat(temp{:,2});
                            end
                        end
                        
                        tempbeh = table2cell(dataDS);
                        
                        numfixes = cellfun(@(x) size(x,1),temptable{:,1},'Uniform',false);
                        numfixes = repmat(numfixes,1,size(tempbeh,2));
                        
                        tempbeh = cellfun(@(x,y) repmat(x,y,1),tempbeh,numfixes,'Uniform',false);
                        
                        tempbeh = cell2table(tempbeh,'VariableNames',dataDS.Properties.VariableNames);
                        
                        
                        common_vars = intersect(tempbeh.Properties.VariableNames,temptable.Properties.VariableNames);
                        
                        if ~isempty(common_vars)
                            tempbeh(:,common_vars) = [];
                        end
                        
                        
                        temptable = horzcat(tempbeh,temptable);
                        
                        dataDS = varfun(@(x) vertcat(x{:}),temptable);
                        dataDS.Properties.VariableNames = temptable.Properties.VariableNames;
                        
                        
                        
                        
                    else
                        %default - we get epoched pupil data
                        if p.Results.fix==false
                            [~,temp]=get(obj,'aligned','event',eyevars{v});
                            fname = eyevars{v};
                            
                            if ~isempty(p.Results.downsample)
                                dataDS.(fname)=downsample(temp(:,2:end)',p.Results.downsample)';
                            else
                                dataDS.(fname)=temp(:,2:end);
                            end
                            
                            dataDS.(strcat(fname,'_window_start')) = repmat(obj.data{1}(1).windows.(fname)(1),height(dataDS),1);
                            dataDS.(strcat(fname,'_window_end')) = repmat(obj.data{1}(1).windows.(fname)(2),height(dataDS),1);
                            
                        else
                           %if fix==true, we get epoched fixation data 
                           
                            for r=1:length(p.Results.rois)
                                [~,temp]=get(obj,'aligned','fix',true,'event',eyevars{v},'roi',p.Results.rois{r});
                                fname = [p.Results.rois{r},'_',eyevars{v}];
                                
                                if ~isempty(p.Results.downsample)
                                    dataDS.(fname)=downsample(temp(:,2:end)',p.Results.downsample)';
                                else
                                    dataDS.(fname)=temp(:,2:end);
                                end
                                
                                
                            end
                            
                            dataDS.(strcat(fname,'_window_start')) = repmat(obj.data{1}(1).fix_windows.(eyevars{v})(1),height(dataDS),1);
                            dataDS.(strcat(fname,'_window_end')) = repmat(obj.data{1}(1).fix_windows.(eyevars{v})(2),height(dataDS),1);
                            
                        end
                        
                        
                end
                
                
                
            end               
        end
        
        
        function obj=remove_artifacts(obj,varargin)
            
            p = inputParser;
            p.addParameter('EL_blinks',1,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('filter_blinks',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('sliding_window',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('blink_window',50,@isnumeric);
            p.addParameter('blink_level',2.0,@isnumeric);
            p.addParameter('outliers',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('outlier_level',2.0,@isnumeric);
            p.addParameter('smooth',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('low_cutoff',1,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('low_cutoff_level',100,@isnumeric);
            p.addParameter('contiguous',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('contiguous_cutoff',2.0,@isnumeric);
            p.addParameter('badxy',1,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('interpolate',0,@(x) islogical(x) || (ismember(x,[0,1])));
            parse(p,varargin{:});
            
            
            numsubs=size(obj.data(:,1),1);
           
            for s=1:numsubs
                 eyeStruct=obj.data{s};
                                            
            
            fprintf('removing artifacts for subject %d out of %d\n',s,numsubs);

            if p.Results.filter_blinks==1 %this needs to be first!
                fprintf('...filtering blinks...\n');
                
                [eyeStruct,~,~] = pupil_filter_blinks(eyeStruct,0);
                
            end
            
            
            if p.Results.EL_blinks==1
                fprintf('...removing blinks based on EyeLink data...\n');
                
                %first remove blinks based on what EyeLink records..
                [eyeStruct,~] = pupil_remove_blinks(eyeStruct);
                
                
                %it's not perfect, so we use some secret sauce..
                %a peak-to-peak sliding window algorithm
                blink_window = p.Results.blink_window * (eyeStruct(1).sample_rate/1000); %window_size is in ms, not samples
                
                if p.Results.sliding_window==1
                    fprintf('...removing artifacts using sliding window (window: %dms, threshold: %3.2f SDs)...\n',p.Results.blink_window,p.Results.blink_level);
                    [eyeStruct, ~] = pupil_sliding_window(eyeStruct,p.Results.blink_window,p.Results.blink_level); %blink_level is in standard deviation units
                
                end
                
                

                
            end
            
 
   
            if p.Results.low_cutoff==1
                fprintf('...removing all samples below %d...\n',p.Results.low_cutoff_level);
            end
            
             if p.Results.outliers==1
                fprintf('...removing outliers (samples beyond %d SDs from the mean)...\n',p.Results.outlier_level);              
                %remove samples that are greater than XX standard
                %deviations from the mean
                if p.Results.outliers==1
                    pupildata=pupil_remove_outliers(eyeStruct,p.Results.outlier_level);
                end
                
            end
             
            if p.Results.badxy==1
                fprintf('...removing bad samples (x,y coords > 10,000)...\n');
            end
            
            
            if p.Results.contiguous==1
                contiguous_cutoff =  p.Results.contiguous_cutoff * (eyeStruct(1).sample_rate/1000); %contiguous_cutoff is in ms, not samples
                fprintf('...removing samples that are not contiguous for at least %dms ...\n',contiguous_cutoff);
            end
            
        
%             padata={eyeStruct.pa}';
%             
%             %             padata=cellfun(@(x) x(eyetracked,:),padata,'UniformOutput',false);
%             maxsamples=max(cellfun(@length,padata));
%             
%             padata=nan(length(eyeStruct),maxsamples);
            
            
            %some functions performed trial by trial here, instead of
            %through helper functions
            for i=1:length(eyeStruct)
                
                pupildata=eyeStruct(i).pa;
               
                %sometimes weird things happen and Eyelink records x/y
                %coordinates as 1,000,000. Remove these samples
                if p.Results.badxy==1
                    
                    badxy = eyeStruct(i).gx>10000 | eyeStruct(i).gy>10000;
                    pupildata(badxy) = NaN;
                    
                    
                end

                %blinks are recorded as zero, but sometimes not exactly
                %zero. Remove anything below 100
                if p.Results.low_cutoff==1
                    pupildata(pupildata<p.Results.low_cutoff_level)=NaN;
                end
                
                
                %for removing little "chunks" that are left after other
                %artifact removals
                if p.Results.contiguous==1
                    %rewrite this with the RunLength function!
                    
                    b=find_contiguous_regions(pupildata); %external function 

                    for j=1:length(b.regwidth)
                        if b.regwidth(j)<p.Results.contiguous_cutoff
                            pupildata(b.regstart(j):b.regend(j))=NaN;
                        end
                        
                    end
                end
                
                
%                 
%                 if length(pupildata)<maxsamples
%                     difference=maxsamples-length(pupildata);
%                     matpupildata=horzcat(pupildata,nan(1,difference));
%                 else
%                     matpupildata=pupildata;
%                 end
%                 
%                 
% %                 padata(i,:)=matpupildata;
                eyeStruct(i).pa=pupildata;
                 
            end
            
            %we interpolate last, to fill in all the gaps
            %by default it will just fill in the gaps. If you choose to
            %smooth, then whole trial is replaced by smoothed data
            if p.Results.interpolate==1 || p.Results.smooth==1
                fprintf('...interpolating using spline fitting...\n');
                if p.Results.smooth==1
                    fprintf('...smoothing data...\n');
                end
                
                eyeStruct = pupil_interpolate(eyeStruct,p.Results.smooth); %2nd argument is for smoothing
                
            end
                
            fprintf('\n\n');
            obj.data{s}=eyeStruct; %save back into the structure
            
            
            end  
            
        end
        
        function obj=filter_trials(obj,varargin)
            p = inputParser;
            p.addParameter('variable','',@ischar);
            p.addParameter('value','');
            p.addParameter('condition','==',@(x) ismember(x,{'==','~=','>=','<=','>','<'}));
            parse(p,varargin{:});
            
            
            
            
            behdata = subsref(obj,struct('type','.','subs',p.Results.variable)); %this sometimes returns cells-within-cells
            
            if iscell(behdata{1})
                
                behdata = cellfun(@(x) vertcat(x{:}),behdata,'Uniform',false); %this pulls everything out of nested cells
            end
           
            keep_idx= cellfun(@(x) eval(['x',p.Results.condition,num2str(p.Results.value)]),behdata,'Uniform',false);
            
            for s=1:length(obj.subs)
                
               obj.data{s} = obj.data{s}(keep_idx{s});
                
            end
            
        end
        
        
        function obj=remove_fields(obj,fieldnames)
            
            
            if ~isa(fieldnames,'cell')
                fieldnames = {fieldnames};
            end
            
            for f=1:length(fieldnames)
                
                for s=1:length(obj.subs)
                    
                    if isfield(obj.data{s},fieldnames{f})
                        obj.data{s} = rmfield(obj.data{s},fieldnames{f});
                    end
                    
                end
                
            end
            
        end
        
        
        function obj=add_fields(obj,varargin)
            p = inputParser;
            p.addParameter('events',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('beh',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('search',{},@iscell);
            p.addParameter('fieldnames',{},@iscell);
            p.addParameter('striptext',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('time',0,@(x) islogical(x) || (ismember(x,[0,1])));
            parse(p,varargin{:});
            
            if p.Results.events
                for s = 1:length(p.Results.search)
                    obj = extract_event(obj,'search',p.Results.search{s},'fieldname',p.Results.fieldnames{s},'striptext',p.Results.striptext,'time',p.Results.time);
                end
                
            elseif p.Results.beh
                for f = 1:length(p.Results.fieldnames)
                    
                     if isfield(obj.data{1},p.Results.fieldnames{f})
                        %if field already exists, remove it first to be
                        %safe
                        obj = remove_fields(obj,p.Results.fieldnames{f});
                    end
                    
                    
                    temp = subsref(obj,struct('type','.','subs',p.Results.fieldnames{f}));
                    
                    for s=1:length(temp)
                        tempcell = num2cell(temp{s});
  
                        [obj.data{s}.(p.Results.fieldnames{f})] = deal(tempcell{:});
                        
                    end
                    
                end
                
                
            end
            
            
        end
        
        
        function obj=extract_event(obj,varargin)
            p = inputParser;
            p.addParameter('search',1,@ischar);
            p.addParameter('fieldname','',@ischar);
            p.addParameter('time',0,@(x) islogical(x) || (ismember(x,[0,1])));
            p.addParameter('striptext',0,@(x) islogical(x) || (ismember(x,[0,1])));
            parse(p,varargin{:});
            
            
            numsubs=size(obj.data(:,1),1);
            
            
            %if no fieldname given, just go with strmatch, minus any
            %unwanted characters
            if isempty(p.Results.fieldname)
                fieldname = regexprep(p.Results.search,'[!@#$%^&()?"*+='',./~` ]','');
            else
                fieldname = p.Results.fieldname;
            end
            
            
            for s=1:numsubs
                
                eyeStruct=obj.data{s};
                
                
                for i=1:length(eyeStruct)
                    
                  
                    find_idx = ~cellfun(@isempty,regexp(eyeStruct(i).events.message,p.Results.search));
                    
                    
                    if ~isempty(find(find_idx,1))
                        
                        if p.Results.time==1
                            val=eyeStruct(i).events.time(find_idx);
                        else
                            val=eyeStruct(i).events.message;
                            val=val{find_idx};
                        end
                        
                        eyeStruct(i).(fieldname) = val;
                    else
                        
                        eyeStruct(i).(fieldname) = '!!MESSAGE_NOT_FOUND!!';
                        
                    end
                    
                    
                end
                
                
                if p.Results.striptext
                    numvector=str2double(regexprep({eyeStruct.(fieldname)}','[a-z A-Z]',''));
                    numvector=num2cell(numvector);
                    [eyeStruct.(fieldname)] = deal(numvector{:});
                end

                
                obj.data{s} = eyeStruct;
           
            end
            
           
            
        end
                       
        
        function obj=rename_messages(obj,search_exp,replace_exp,replacement)
           
            numsubs=size(obj.data(:,1),1);
           
            for s=1:numsubs
                
                eyeStruct=obj.data{s};
            
                for i = 1:length(eyeStruct)
                   
                    messages = eyeStruct(i).events.message;
                    matches =  ~cellfun(@isempty,regexp(messages,search_exp));
                    
                    if ~isempty(find(matches,1))
                        
                        eyeStruct(i).events.message(matches) = regexprep(messages(matches),replace_exp,replacement);
                    end
                    
                    
                    
                end
                
                
                obj.data{s} = eyeStruct;
                
                
            end
            
            
            
        end
        
        function obj=add_event(obj,times,varargin)
            p = inputParser;
            p.addParameter('index',{},@iscell);
            p.addParameter('name','',@ischar);
            parse(p,varargin{:});
            
            switch class(times)
                
                case {'table'}
                    
                    times = times(:,horzcat(obj.subject_var,p.Results.index,p.Results.name));
                
                case {'double','single'}
                    times = array2table(times,'VariableNames',{obj.subject_var,p.Results.name});                    
                     
                case {'cell'}
                    times = vertcat(times{:});
                    times = cell2table(times,'VariableNames',{obj.subject_var,p.Results.name});
                    
            end
            
            
           %the work is done with add_behdata because the process is the
           %same. 
           obj = add_behdata(obj,times,'index',p.Results.index,'addevent',true,'append',true);
            
            
            
        end
        
        function obj=load_events(obj,varargin)
            
            
            numsubs = length(obj.subs);
            
            for s=1:numsubs
                
                eyeStruct = obj.data{s};
                
                for t =1:length(eyeStruct)
                
                    eyeStruct(t).events=getEdfMessages(eyeStruct(t));%external file
            
                end
                
                obj.data{s} = eyeStruct;
                
                
            end
            
        end
      
        
        function obj=pupil_subset(obj,criteria)
           
            numsubs=size(obj.data(:,1),1);
            
            expr=sprintf('beh.%s & ',criteria{:});
            expr=expr(1:end-2);
            
            for s=1:numsubs
               
                beh=obj.data{s,2};
                eyes=obj.data{s};
                
                eval(['keep=' expr ';']);
                
                beh=beh(keep,:);
                eyes=eyes(keep);
                
                
                obj.data{s,2}=beh;
                obj.data{s}=eyes;
                
            end
            
        end
      
        function [fixcell,fixtable] = get_Fixations(obj,varargin)
           
            p = inputParser;
            p.addParameter('subjects',[],@(x) isnumeric(x) || iscell(x));
            p.addParameter('fields',{},@(x) iscell(x))
            parse(p,varargin{:});
            
            
            %get all fixation data at once
            S = struct;
            S.type='.';
            S.subs='Fixations';
            
            %this gives a cell for each subject, with a cell array for each
            %trial 
            temp = subsref(obj,S);
            
            fnames = fieldnames(temp{1}{1});
            
            %convert the nested cell arrays to matricies
            fixcell = cellfun(@(y) cellfun(@struct2array,y,'Uniform',false),temp,'Uniform',false);
            
            
            %if we only want certain fields, this pulls out the right ones
            %(in the order we want)
            if ~isempty(p.Results.fields)
                f_idx = cell2mat(cellfun(@(x) find(strcmp(x,fnames)),p.Results.fields,'Uniform',false));
                fixcell = cellfun(@(y) cellfun(@(z) z(:,f_idx),y,'Uniform',false), fixcell,'Uniform',false);
                fnames = p.Results.fields;
            end
            
            
            %now to produce the table, we have to extract the subject and
            %trial data
            
            %this gives an array of 1:numtrials for each subject
            trialnums = cellfun(@(y) num2cell(transpose(1:length(y))),fixcell,'Uniform',false);
            
            %this replicates the trial number by the number of fixations
            %for each trial
            trialnums = cellfun(@(x,y) cellfun(@(z,q) repmat(q,size(z,1),1),x,y,'Uniform',false),fixcell,trialnums,'Uniform',false);
            
            %this makes it into a vector (instead of cells within cells)
            trialnums = cellfun(@(x) vertcat(x{:}),trialnums,'Uniform',false);

            %this replicates the subject ID for the number of fixations for
            %each subject
            subids = cellfun(@(x,y) repmat(y,size(x,1),1),trialnums,obj.subs,'Uniform',false);
            
            %this makes everything into a big vector
            subids = vertcat(subids{:});
            trialnums = vertcat(trialnums{:});
            
            
            %stacks fixation data as a big matrix
            temp = cellfun(@(y) vertcat(y{:}),fixcell,'Uniform',false);
            temp = vertcat(temp{:});
            
            %add the subject and trial information, then convert to a table
            temp = horzcat(subids,trialnums,temp);
            fixtable = array2table(temp,'VariableNames',horzcat(obj.subject_var,'trial_idx',fnames));
                    
        end
        
     
        function obj=index_fixations(obj)
            
           %grab all fixations 
           fixcell = get_Fixations(obj,'fields',{'sttime','entime'});
           
           %for each subject
           for s = 1:length(fixcell)
              
               fixdata =fixcell{s,1}; %pull out the data
               
               for trial =1:size(fixdata,1) %loop through each trial
                    %this is all zeros at first (#fixations x #timepoints)
                    idxmat = obj.data{s}(trial).fixation_times;
                    
                    for fix=1:size(idxmat,1) %loop through each fixation
                       
                      if ~isnan(fixdata{trial}(fix,1)) 
                       idxmat(fix,fixdata{trial}(fix,1):fixdata{trial}(fix,2)) = 1; %mark the timepoints where the fixations occurred
                      end
                        
                    end
                   
                    obj.data{s}(trial).fixation_times = idxmat; %save it
                    
                    fix_mask = max(idxmat,[],1);
                    fix_mask(fix_mask==0) = NaN;
                    obj.data{s}(trial).fixation_time_mask = fix_mask;
                   
               end
                    
           end
            
            
            
        end
      
        
        function obj=clearROIs(obj)
           
            obj.rois.single=[];
            obj.rois.combined=[];
            
        end
        
        
        function obj=makeROIs(obj,pos,varargin)
            
            p = inputParser;
            p.addParameter('radius',50,@isnumeric);
            p.addParameter('shape','circle',@(x) ismember(x,{'circle','circular','ellipse','elliptical','square'}));
            p.addParameter('xradius',50,@isnumeric);
            p.addParameter('yradius',10,@isnumeric);
            p.addParameter('angle',0,@(x) min(x)>=0 && max(x)<= 360);
            p.addParameter('clear',0);
            p.addParameter('names',{},@iscell);
            parse(p,varargin{:});
            

            if p.Results.clear==1
                obj = clearROIs(obj);
            end
            
            xres = obj.screen.dims(1);
            yres = obj.screen.dims(2);
            

            if isempty(obj.rois.single)
                existing_rois = 0;
            else
                existing_rois = length(obj.rois.single);
            end
            
            
            if length(p.Results.angle)==1
                angles = repmat(p.Results.angle,size(pos,1),1);
            else
                angles=p.Results.angle;
            end
            
            
            if isempty(p.Results.names)
%                 names = num2cell((existing_rois+1):(existing_rois+size(pos,1)));
                temp = (existing_rois+1):(existing_rois+size(pos,1));
                names = strread(num2str(temp),'%s');
            else
                names = p.Results.names;
            end
            
            
            for r =(existing_rois+1):(existing_rois+size(pos,1))
            
             i = r - existing_rois;  %index of current set of rois 
             xpos = pos(i,1);
             ypos = pos(i,2);
             
             xcenter = fix(obj.screen.dims(1)/2);
             ycenter = fix(obj.screen.dims(2)/2);
             
             [XX, YY] = meshgrid(0:(xres-1),0:(yres-1));
             
             obj.rois.single(r).name = names{i};
             obj.rois.single(r).coords = [xpos,ypos];
             obj.rois.single(r).shape = p.Results.shape;
             
             switch p.Results.shape
                 case {'circle','circular'}
                     obj.rois.single(r).mask = sqrt((XX-xpos).^2+(YY-ypos).^2)<=p.Results.radius;
                     
                 case {'ellipse','elliptical'}
                    
                    if angles(i)>0
                     
                     xshift = xpos - xcenter;
                     yshift = ypos - ycenter;
                        
                     %create an ellipse in the center, then rotate   
                     el=((XX-xcenter)/p.Results.xradius).^2+((YY-ycenter)/p.Results.yradius).^2<=1;   
                     el=imrotate(el,angles(i),'nearest','crop');
                     
                     %then shift the image so it's centered over the
                     %correct point
                     RA = imref2d(size(el)); %so we keep the same image size
                     tform = affine2d([1 0 0; 0 1 0; xshift yshift 1]);
                     el = imwarp(el, tform,'OutputView',RA);
                     
                    else
                        el=((XX-xpos)/p.Results.xradius).^2+((YY-ypos)/p.Results.yradius).^2<=1;
                    end
                    
                     obj.rois.single(r).mask =el;
             end
             
            end
            
        end
        
        function obj=combineROIs(obj)
            
            numrois = length(obj.rois.single);
            
            combined = zeros(size(obj.rois.single(1).mask));
            
            for r = 1:numrois
               
                combined = combined + obj.rois.single(r).mask;  
 
            end
            
            obj.rois.combined = combined;
            
        end
        
       
        function obj= calcFixationHits(obj,varargin)
            
            p = inputParser;
            p.addParameter('rois','all',@(x) iscell(x) || ischar(x));
            parse(p,varargin{:});
            
              
            numsubs = length(obj.subs);
            
            
            if ~iscell(p.Results.rois) && strcmp(p.Results.rois,'all')
                rois = {obj.rois.single.name};
            elseif ~iscell(p.Results.rois) 
                rois = {p.Results.rois};
            else
                rois = p.Results.rois;
            end
  
            numROIs = length(rois);
            
            allfixdata = get_Fixations(obj,'fields',{'gavx','gavy'});
            
            for s=1:numsubs
                
                fixation_hits = struct;
                data = allfixdata{s};
                
                
                for r=1:numROIs
           
                    xres = obj.screen.dims(1);
                    yres = obj.screen.dims(2);
                    
                    roi_idx = find(ismember({obj.rois.single.name},rois{r}));
                    
                    roi_mask = obj.rois.single(roi_idx).mask;
                    
                

                    
                    %%
                    %matrix of all coordinates
                    coords = double(vertcat(data{:}));
                    
                    
                    %creates a vector of the trial #, repeated for each
                    %fixation in that trial. Useful for putting our data
                    %back into the cell array format
                    trialnums = num2cell(1:length(data))';
                    numfixes = cellfun(@(x) size(x,1),data(:,1),'Uniform',false);
         
                    trial_vec = cellfun(@(x,y) repmat(y,x,1),numfixes,trialnums,'Uniform',false);    
                    trial_vec = vertcat(trial_vec{:});

                    %bad samples recoded as zero
                    coords(coords(:,1)>xres | coords(:,1)<=0 | isnan(coords(:,1)),1) = xres;
                    coords(coords(:,2)>yres | coords(:,2)<=0 | isnan(coords(:,2)),2) = yres;

                    coords = ceil(coords);
                    
                    %find all unique fixations
                    [ucoords,idx_u,idx_a] = unique(coords,'rows','stable');

                    scr = zeros(yres,xres);

                    %find the indicies of the fixations
                    idx = sub2ind([768,1024],ucoords(:,2),ucoords(:,1));
                    scr(idx)=1:length(ucoords); 

                    fix_bin = logical(scr); %binary mask of fixations
                    overlap = (double(roi_mask+fix_bin)>1); %the overlap between the mask and the fixations

                    overlap = double(overlap).*scr;

                    overlap_idx = unique(overlap);
                    overlap_idx = overlap_idx(overlap_idx>0); %rows in ucoords that overlap with mask
                    
                    %now map this back to the original data (all fixations,
                    %not just the unique ones)
                    hits = idx_a(ismember(idx_a,overlap_idx));
                    hits_sub =zeros(length(ucoords),1);
                    hits_sub(hits)=1;
                    
                    %binary vector-- hit or not
                    hits_all=single(hits_sub(idx_a));

                    if isnumeric(rois{r})
                        name = strcat('roi_',num2str(rois{r}));
                    else
                        name = strcat('roi_',obj.rois.single(roi_idx).name);
                    end


                    %take our 1 big vector and split to a cell to put back
                    %in our array (1 cell for each trial)
                    hits_all = split_by_idx(hits_all,trial_vec);


                    [fixation_hits(1:length(hits_all)).(name)] = deal(hits_all{:});
                    
                end
                
                %convert our fixation_hits structure to a cell. Why?
                %because matlab is dumb.
                temp = arrayfun(@(x) {x},fixation_hits);
                
                %now we can fill in our data
                [obj.data{s}.fixation_hits] = deal(temp{:}); 
                
            end
               
        end
        
        
        function varargout = subsref(obj,S)
            %this allows you to index the object like you would with
            %structures. You can use any field in obj.data or any behavioral
            %variable in obj.data{i}.beh
            
            
            if any(strcmp(S(1).type,'.')) && ismember(S(1).subs,fieldnames(obj));
                %if we're dealing with the highest-level fields, just use
                %matalab's builtin stuff
                [varargout{1:nargout}] = builtin('subsref',obj,S);
            
            elseif length(S)==1
                [objfields,behfields] = get_all_fields(obj);
                
                if ismember(S.subs,objfields)
                   [varargout{1:nargout}] = cellfun(@(x) reshape({x.(S.subs)},[],1),obj.data,'Uniform',false);
%                    [varargout{1:nargout}] = cellfun(@(x) vertcat(x{:}),allobj,'Uniform',false);
                   
                elseif ismember(S.subs,behfields)
                     
                    allbeh = cellfun(@(x) [x.beh],obj.data,'Uniform',false);
                    [varargout{1:nargout}] = cellfun(@(x) reshape([x.(S.subs)],[],1),allbeh,'Uniform',false); 
                else
                    error('field "%s" not found!',S(1).subs)
                end
            else
                error('Sorry, cannot do multi-level references (yet)');
                
            end
            
        end
        
        
        function [objfields,behfields] = get_all_fields(obj)
            
            objfields = cellfun(@fieldnames,obj.data,'Uniform',false);
            
            if isfield(obj.data{1},'beh')
                behfields =  cellfun(@(x) fieldnames(x(1).beh),obj.data,'Uniform',false);
            else
                behfields = {''};
            end
                 
            objfields = vertcat(objfields{:});
            behfields = vertcat(behfields{:});
            
            objfields = reshape(objfields,[],1);
            behfields = reshape(behfields,[],1);
            
            objfields = unique(objfields);
            behfields = unique(behfields);
             
        end
        
        
        
        function obj=mapROIs(obj,newName,varargin)
            p = inputParser;
            p.addParameter('rois','all',@iscell);
            p.addParameter('indicator',{});
            parse(p,varargin{:}); 
            
            numsubs = length(obj.subs);
            
            %indicator should be a cell array, 1 cell per subject
            %each cell has a row for each trial, with either a number or
            %the name of the roi
            
            
            if ischar(p.Results.indicator)
                %if variable name given, pull that information out
                idx = subsref(obj,struct('type','.','subs',p.Results.indicator));
            else
                %otherwise, take the input as given
                idx = p.Results.indicator;
            end
            
            %if it's the same for each subject, can give 1 cell array.
            %Repeat it for each subject
            if length(idx)==1 && length(obj.subs)>1
                idx = repmat(idx,length(obj.subs),1);
            end
            
            by_name = false;
            
            %if rois specified by name (otherwise, assume they're specified
            %by index)
            if isa(idx{1},'cell')
                by_name = true;
            end
            
            
            newName = strcat('roi_',newName);
            
            for s=1:numsubs

               
               roimap = idx{s};

               for i=1:length(obj.data{s})
                   %loop through each trial, get the name of the original
                   %rois that match up with the trial-specific one, then
                   %just copy the data
                   
                   
                   if by_name %if specified by name
                        roiname = strcat('roi_',roimap{i});
                   else
                       %find the corresponding roi based on the index
                       %number
                       if roimap(i) > 0 && roimap(i) <= length(obj.rois.single)
                            roiname = strcat('roi_',obj.rois.single(roimap(i)).name);
                       else
                           roiname = [];
                       end
                   
                   if ~isempty(roiname)
                        obj.data{s}(i).fixation_hits.(newName) = obj.data{s}(i).fixation_hits.(roiname);
                   else
                       obj.data{s}(i).fixation_hits.(newName) = single(nan(size(obj.data{s}(i).fixation_times,1),1));
                   end
                   
                   
                   end
               end
                
                
            end
            
            
            
        end
        
        
        function obj=calcHitsOverTime(obj,varargin)
            p = inputParser;
            p.addParameter('rois','all',@iscell);
            parse(p,varargin{:});
            
            numsubs = length(obj.subs);
            
            
            if ~iscell(p.Results.rois) && strcmp(p.Results.rois,'all')
                rois = {obj.rois.single.name};
            elseif ~iscell(p.Results.rois)
                rois = {p.Results.rois};
            else
                rois = p.Results.rois;
            end
            
            numROIs = length(rois);
            
            fixation_hits_all =  subsref(obj,struct('type','.','subs','fixation_hits'));
            fixation_times_all = subsref(obj,struct('type','.','subs','fixation_times'));
            
            fixation_time_mask_all = subsref(obj,struct('type','.','subs','fixation_time_mask'));
            
            for s=1:numsubs
                
                fixation_hits = fixation_hits_all{s};
                fixation_times = fixation_times_all{s};
                fixation_time_mask = fixation_time_mask_all{s}; %shows where we have no fixation data

                [obj.data{s}.fix] = deal([]);
                
                allhits = struct;
                
                for r=1:numROIs
                    

                    roihits = cellfun(@(x) x.(strcat('roi_',rois{r})),fixation_hits,'Uniform',false);
                    
                    hitsovertime = cellfun(@(x,y) double(x)'*double(y),roihits,fixation_times,'Uniform',false);
                    
                    %multiply by the mask so that missing data doesn't
                    %contribute to the average (times where no fixation
                    %will be NaN instead of zero)
                    hitsovertime = cellfun(@(x,y) double(x) .* double(y), hitsovertime,fixation_time_mask,'Uniform',false);
                    
                    [allhits(1:length(hitsovertime)).(rois{r})] = deal(hitsovertime{:}); 

                    
                end
                
                %convert our allhits structure to a cell. Why?
                %because matlab is dumb
                temp = arrayfun(@(x) {x},allhits);
                
                [obj.data{s}.fix] = deal(temp{:});

                
            end
            
            
        end
        
        
            
        function aggDS=aggregate(obj,varargin)
            p = inputParser;
            p.addParameter('fix',true);
            p.addParameter('factors',{},@iscell);
            p.addParameter('eyedata',{},@iscell);
            p.addParameter('rois',{},@iscell);
            p.addParameter('function',@nanmean);
            p.addParameter('filter','');
            p.addParameter('value',1);
            parse(p,varargin{:});
            
            
            if ~isempty(p.Results.filter)
                dataDS=build_dataset(obj,horzcat(obj.subject_var,p.Results.factors,p.Results.filter),p.Results.eyedata,'fix',p.Results.fix,'rois',p.Results.rois); %build behavior-pupil dataset
                dataDS= dataDS(dataDS.(p.Results.filter)==p.Results.value,:);
            else
            
                dataDS=build_dataset(obj,horzcat(obj.subject_var,p.Results.factors),p.Results.eyedata,'fix',p.Results.fix,'rois',p.Results.rois); %build behavior-pupil dataset
            end
            


            
            dataVarnames = dataDS.Properties.VariableNames((end-length(p.Results.rois)+1):end);
            
            %first aggregate by subject
            aggDS = grpstats(dataDS,horzcat(obj.subject_var,p.Results.factors),p.Results.function,'DataVars',dataVarnames);
            aggDS.Properties.VariableNames((end-(length(p.Results.rois))+1):end) = dataVarnames;
            
            %then across subjects
            aggDS = grpstats(aggDS,p.Results.factors,p.Results.function,'DataVars',dataVarnames);
            aggDS.Properties.VariableNames((end-(length(p.Results.rois))+1):end) = dataVarnames;
            
             
        end
        
        function scatterplots(obj,separate_by,varargin)
            p = inputParser;
            p.addParameter('rois',{},@iscell);
            p.addParameter('overlay',false);
            parse(p,varargin{:});
           
           alldata = build_dataset(obj,'beh',separate_by,'allfix',true,'rois',horzcat('xy',p.Results.rois));
           alldata= alldata(~any(ismissing(alldata),2),:); %remove rows with any missing data
           
           
           if ~isempty(p.Results.rois)
               separate_by = horzcat(separate_by,strcat(p.Results.rois,'_hit'));
           end
           
           
           
           fixScatterplots(obj,alldata,separate_by,'x_coord','y_coord');
           
          placeFigures('intergapH',10,'intergapV',10,'toolsize',30)
            
            
        end
         
        
        function fixScatterplots(obj,fixdata,separate_by,x_col,y_col,varargin)
            %Specify a dataset, factors to separate the data by (in a cell array), and x and y
            %columns (as variable names), and it will split into the appropriate number of plots. It's a
            %recursive function that handles an arbitrary number of factors. Factors
            %should be listed "hierarchically"-- the last one listed will group into
            %subplots on a single figure. Otherwise, they will be separated into
            %separate figures.
            %
            %
            %Example: if we have 5 subjects and 2 conditions each,
            %fixationScatter(testDS,{'SUBID' 'CONDITION'},'X','Y') generates 5 figures with 2 plots each.
            %fixationScatter(testDS,{'CONDITION' 'SUBID'},'X','Y') generates 2 figures
            %with 5 subplots each.
            
            
            
            
            
            %Title is optional when fucntion is called (but fed into further recursive
            %steps). Initialize as blank.
            bigtitle='';
            
            
            %initialize to this, if bigtitle is provided, it will be used as the
            %filename
            filename=[separate_by{1} '.tif'];
            
            %check for optional x and y limits for plots e.g., xlimits=[0 1024]
            if (length(varargin)>=2 && ~isempty(varargin{1}))
                xlimits=varargin{1};
                ylimits=varargin{2};
            else
                xlimits = [0,obj.screen.dims(1)];
                ylimits = [0,obj.screen.dims(2)];
            end
            
            %check for optional title argument
            if (length(varargin)>=3 && ~isempty(varargin{3}))
                bigtitle=varargin{3};
                filename=[regexprep(bigtitle,'=', '-') '.tif'];
            end
            
            %check for optional overlay argument (3xN matrix, specifying x,y, and
            %radius of each circle)
            if (length(varargin)>=4)
                overlay_coords=varargin{4};
            end
            
            if (length(varargin)>=5)
                movieMode=varargin{5};
            else
                movieMode=0;
            end
            
            %grab the fist factor from the list
            groupvar=separate_by{1};
            
            %base condition, if we only have 1 factor
            if length(separate_by)==1
                
                %figure out how many subplots we need
                all_conditions = unique(fixdata.(groupvar));
                plotrows = ceil(sqrt(length(all_conditions)));
                
                %create new figure
                big_figure=figure();
                                    
                   
                %for each factor level...
                for i=1:length(all_conditions)
                    %subset data...
                    current = all_conditions(i);
                    groupvar=separate_by{1};
                    current_subset = fixdata(fixdata.(groupvar)==current,:);
                    
                    %...create new subplot (2 plots side-by-side if there are only 2, otherwise
                    %make a square matrix of plots (e.g., 3x3 for 9 or 10 plots).
                    
                    
                    if(~movieMode && length(all_conditions)==2)
                        subplot(1,2,i)

                    
                    elseif (~movieMode && length(all_conditions)~=2)
                        subplot(plotrows,plotrows,i);

                    end
                    
                    if movieMode
                        set(gca,'nextplot','replacechildren');
                        xlim(obj.screen.dims(1));
                        ylim(obj.screen.dims(2));
                    end
                    
                    %draw the scatterplot of the current subset (red dots)
                    h = scatter(current_subset.(x_col),current_subset.(y_col),10,'r');
                    hold on;
                    
                    
                    %if we need to draw overlays, draw them, using helper drawOverlay
                    %function
                    if (exist('overlay_coords'))
                        for i=1:size(overlay_coords,1)
                            drawOverlay(overlay_coords(i,:));
                            hold on;
                        end
                    end
                    
                        set(gcf,'Toolbar','none','Color','white','Menubar','none');
    
                        set(gca,'LooseInset',get(gca,'TightInset'));
                        
                        set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
                    
                    %title for each subplot, written below (so it doesn't overlap with
                    %overall title.
                    title_text = sprintf('%s = %d',groupvar,current);
                    title_text = strrep(title_text,'_','-');
                    xlabel(title_text);
                    
                    %set x and y limits, if given
                    if (exist('xlimits') && exist('ylimits'))
                        ylim(ylimits);
                        xlim(xlimits);
                    end
                    
                    if ~movieMode
                        %create 1 big title at the top of the figure
                        figure(big_figure);
                        set(big_figure,'NextPlot','add');
                        axes;
                        h = title(strrep(bigtitle,'_','-'));
                        set(gca,'Visible','off');
                        set(h,'Visible','on');
                    end
                    

                    
                    
                    if movieMode
                        mov=getframe(big_figure);
                        pause(.1);
                        imwrite(mov.cdata,filename,'writemode','append');
                    end
                    
                end
                
                
                %otherwise, do a recursive function call...
            else
                
                %we have to pass the optional arguments this time, so just set them as
                %empty if they weren't given the first time.
                if ~exist('xlimits')
                    xlimits=[];
                end
                
                if ~exist('ylimits')
                    ylimits=[];
                end
                
                if ~exist('overlay_coords')
                    overlay_coords=[];
                end
                
                %for each level of the first factor...
                all_conditions = unique(fixdata.(groupvar));
                for i=1:length(all_conditions)
                    %...create a subset
                    current = all_conditions(i);
                    current_subset = fixdata(fixdata.(groupvar)==current,:);
                    
                    %...give it an informative title.
                    title_text = sprintf('%s %s = %d',bigtitle,groupvar,current);
                    
                    %...and call fixationScatter again, this time with the subset and
                    %with the remaining factors.
                    
                    fixScatterplots(obj,current_subset,separate_by(2:end),x_col,y_col,xlimits,ylimits,title_text,overlay_coords,movieMode);
                    
                end
            end
            
        end
       
        function obj = binData(obj,window_size,varargin)
            
            p = inputParser;
            p.addParameter('events',{},@iscell);
            p.addParameter('fix',true);
            p.addParameter('rois',{'none'},@iscell);
            p.addParameter('func','mean',@ischar);
            parse(p,varargin{:});
            
            
            if obj.data{1}(1).sample_rate ~= 1000
                %window size is in ms, this makes sure that stays true
                window_size = window_size*(z.data{1}(1).sample_rate/1000);
                
            end
            
            for e = 1:length(p.Results.events)
                %loop through each event
                
                for r = 1:length(p.Results.rois)
                    %and for each roi if it's fixation data
                    
                    if p.Results.fix
                        aligned = get(obj,'aligned','fix',p.Results.fix,'event',p.Results.events{e},'roi',p.Results.rois{r});
                        aligned = aligned(:,2); %the first column is just the subject ids.
                        
                    else
                        aligned = get(obj,'aligned','event',p.Results.events{e});
                        aligned = aligned(:,2);
                        
                    end
                    
                    binned = cell(size(aligned));
                    
                    for s = 1:length(aligned)
                        %for each subject...
                        
                        temp = aligned{s};
                        
                        if p.Results.fix
                            temp(isnan(temp)) = -Inf; %so our running maximum algorithm works..
                            
                            temp2 = downsample(minmaxfilt(temp,[1,window_size],'max','valid')',window_size)'; %minmaxfilt is from the file exchange, does a running maximum (quickly).
                            %downsample pulls out every x samples (which corresponds to our
                            %bins).
                            
                            temp2(isinf(temp2)) = NaN; %put the NaNs back
                        else
                            %if pupil data, then do a moving average or moving
                            %median
                            switch p.Results.func
                                case {'mean'} %default
                                    temp2 = downsample(moving_average(temp,window_size,'valid')',window_size)';
                                case {'median','med'}
                                    temp2 = downsample(medfilt2(temp,[1,window_size])',window_size)';
                            end
                        end
                        
                        
                        binned{s}= temp2;
                        
                    end
                    
                    %put the data into the object to replace the original
                    obj = put(obj,binned,'fix',p.Results.fix,'event',p.Results.events{e},'roi',p.Results.rois{r});
                    
                    
                end
                
            end
        end
       
        
        function obj = put(obj,input,varargin)
            p = inputParser;
            p.addParameter('event','',@ischar);
            p.addParameter('fix',true);
            p.addParameter('roi','',@ischar);
            parse(p,varargin{:});
            
            
            if length(input)~=length(obj.subs)
                error('size of input data does not match number of subjects!');
            end
            
            
            for s=1:length(obj.subs)
               
                if p.Results.fix
                    
                    for t = 1:length(obj.data{s})
                       
                        obj.data{s}(t).fix_aligned.(p.Results.event).(p.Results.roi) = input{s}(t,:);
                        
                        
                    end
                    
                else
                    for t = 1:length(obj.data{s})
                        
                        obj.data{s}(t).aligned.(p.Results.event) = input{s}(t,:);
                        
                        
                    end
                
                end
                
                
            end
            
            
        end
        
        
    end
    
end
    
    
   









