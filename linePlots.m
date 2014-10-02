function allFigData=linePlots(allData,varargin)

p = inputParser;
p.addParameter('plotVars',{},@iscell);
p.addParameter('lineVars',{},@iscell);
p.addParameter('dataVar','',@ischar);
p.addParameter('subVar','ID',@ischar);
p.addParameter('fix',true,@(x) islogical(x) || ismember(x,[0,1]));
p.addParameter('roi','',@ischar);
p.addParameter('legend',true,@(x) islogical(x) || ismember(x,[0,1]))
p.addParameter('LineWidth',2,@isnumeric)
p.addParameter('ylim','auto',@(x) (length(x)==2 && isnumeric(x)) || strcmp(x,'auto'));
p.addParameter('ttest',false,@(x) islogical(x) || ismember(x,[0,1]))
p.addParameter('alpha',.05,@isnumeric)
p.addParameter('filter','',@ischar);
p.addParameter('labels',{},@iscell);
p.addParameter('func',@nanmean);
p.addParameter('value',1);

parse(p,varargin{:});

lineVars = p.Results.lineVars;
factors = p.Results.plotVars;
dataVar = p.Results.dataVar;

cname = 'Set2';
ctype = 'qual';


if ~isempty(p.Results.filter)
    allData = allData(allData.(p.Results.filter)==p.Results.value,:);
end




if isempty(factors) %if you just want 1 plot with multiple lines
    allData.ONEPLOT = ones(length(allData.(dataVar)),1);
    factors = {'ONEPLOT'};
end

[allFigData,group_idx,group_levels]=makeSubsets(allData,factors);



for i=1:length(allFigData)
    
    allFigData(i).subset.lineVars=makeLabels(allFigData(i).subset,lineVars);
    
    %aggregate within subject first
    temp=grpstats(allFigData(i).subset,{p.Results.subVar,'lineVars'},p.Results.func,'DataVars',dataVar);
    temp.Properties.VariableNames{end}=dataVar;
    
    if ~isempty(lineVars) %if we need to aggregate data
        
        temp=grpstats(temp,'lineVars',@nanmean,'DataVars',dataVar);
        temp.Properties.VariableNames{end}=dataVar;
        
      
        
    else %otherwise, plot all the lines
        
        allFigData(i).subset.lineVars=num2str([1:height(allFigData(i).subset)]');
        temp=allFigData(i).subset;
        
    end
    
    temp = sortrows(temp,'lineVars');
    
    allFigData(i).subVar=p.Results.subVar;
    allFigData(i).dataVar=dataVar;
    
    allFigData(i).plotDS=temp;
    allFigData(i).plotdata=double(temp.(dataVar));
    
    
    h=plotRows(allFigData(i).plotdata,ctype,cname,p.Results.LineWidth,p.Results.ylim);
    
    
    
    allFigData(i).FIG_H=h;
    allFigData(i).AH=findobj(h,'Type','axes');
    
    if p.Results.legend
        names=temp.lineVars;
        clickableLegend(cellfun(@fix_title_string,names,'Uniform',false),'Location','NorthWest');
    end
    
    title(cellfun(@fix_title_string,allFigData(i).title,'Uniform',false),'fontSize',10);
    
    set(h,'Name',allFigData(i).title{:},'Toolbar','none','Color','white','Menubar','none');
    
    set(allFigData(i).AH,'LooseInset',get(gca,'TightInset'));
    
%     set(findall(h,'type','text'),'fontSize',12)
    
    D=get(gca,'Children'); %get the handle of the line object
    allFigData(i).fig_data=get(D);
    
    allFigData(i).lineVars=lineVars;
          
end


set(findobj('Type','axes'),'box','off');


screensize=get(0,'ScreenSize');
minwidth=floor(screensize(4)/(round(length(allFigData)/2)));
minwidth=max(minwidth,150);

placeFigures('intergapH',10,'intergapV',10,'toolsize',30);

if p.Results.ttest && length(p.Results.lineVars)==1
    allFigData=dottests(allFigData,p.Results.alpha);
end

%useful plot stuff
% set(findobj('Type','line'),'LineWidth',4)

end

%%

function h=plotRows(data,ctype,cname,linewidth,ylimits)

h=figure;

if isa(ylimits,'char') && strcmp(ylimits,'auto')
    ylimits=[min(data(:))-min(data(:))*.20, max(data(:))*1.2];
end


% cmap=colormap(hsv(size(data,1)));
ctype=cbrewer(ctype,cname,max(size(data,1),3));

for i=1:size(data,1)
    y=data(i,:);
    x=1:size(y,2);
    
    plot(x,y,'Linewidth',linewidth,'color',ctype(i,:)); hold all;
    
end


ylim(ylimits);


hold off;

return


end


function allfigdata=dottests(allfigdata,alphalevel)

allData=vertcat(allfigdata(:).subset);

subVar=allfigdata(1).subVar;
lineVars=allfigdata(1).lineVars{1};
dataVar=allfigdata(1).dataVar;

levels=unique(allData.(lineVars));

if length(levels)>2
    error('Cannot perform ttests on factors with more than 2 levels!');
end

for i=1:length(allfigdata)
    
   temp=grpstats(allfigdata(i).subset,{subVar,'lineVars'},@nanmean,'DataVars',dataVar);
   temp.Properties.VariableNames{end}=dataVar; 
   
   temp.lineVars=grp2idx(temp.lineVars);
    
   level1=temp{temp.lineVars==1,dataVar};
   level2=temp{temp.lineVars==2,dataVar};
   
   sig=nan(1,size(temp.(dataVar),2));
   
   for j=1:size(temp.(dataVar),2)
      
       [sig(j),~,~]=ttest(level1(:,j),level2(:,j),'Alpha',alphalevel);
       
   end
   
   allfigdata(i).ttest=sig;
   
   
   plotdata=allfigdata(i).plotdata;
   
   minmax=nan(size(plotdata));
   
   for k=1:size(plotdata,2)
      minmax(1,k)=min(plotdata(:,k));
      minmax(2,k)=max(plotdata(:,k)); 
   end
   
   
   sig(sig==0)=NaN;
   if ~all(isnan(sig))
       sigRegions=find_contiguous_regions(sig);
       
       for s=1:length(sigRegions)
           
           interval = sigRegions.regstart(s):sigRegions.regend(s);
           figure(allfigdata(i).FIG_H);
           fillRegion(interval,minmax(2,interval),minmax(1,interval),rgb('orange'),.2);
       end
       
   end
   
   set(findobj('Type','axes'),'Visible','on');
end


end


function h=fillRegion(xpoints,upper,lower,color,transparency)

if length(upper)==length(lower) && length(lower)==length(xpoints)
    filled=[upper,fliplr(lower)];
    xpoints=[xpoints,fliplr(xpoints)];

    hold on

    
    fillhandle=fill(xpoints,filled,color);%plot the data

    set(fillhandle,'EdgeColor',color,'FaceAlpha',transparency,'EdgeAlpha',0);%set edge color
    

    hold off

    uistack(fillhandle,'bottom');
    set(findobj('Type','axes'),'box','off');

end
end
