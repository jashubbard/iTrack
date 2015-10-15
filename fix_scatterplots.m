function fix_scatterplots(ds,factors,varargin)
%given some fixation data, creates separate scatterplots for each
%combination of different factors. Can also add an overlay (binary mask the
%same size as the screen, with 1 representing a given ROI). This is the
%function called by the iTrack scatterplots and quickview functions


p = inputParser;
p.addParameter('zoom',true,@(x) islogical(x) || ismember(x,[0,1]));
p.addParameter('screen',true,@(x) islogical(x) || ismember(x,[0,1]));
p.addParameter('overlay',[],@(x) islogical(x) | isnumeric(x));
p.addParameter('hitonly',false,@(x) islogical(x) || ismember(x,[0,1]));
p.addParameter('skipnan',true,@(x) islogical(x) || ismember(x,[0,1]));
p.addParameter('size',4,@isnumeric);
p.addParameter('Xvar','fix_x',@ischar);
p.addParameter('Yvar','fix_y',@ischar);
p.addParameter('color','yellow',@ischar);
p.addParameter('roi',{},@iscell);
p.addParameter('screendims',[1024,768],@isnumeric);

parse(p,varargin{:});


%get fixation x/y coordinates
%             fixdata = build_dataset(obj,'allfix',true,'rois',{'xy'});

%             fixdata = get_new(obj,'fixations');

%             fixdata = report(ds,'fixations','behvars',factors,'rois',p.Results.roi);



%create a factor for whether the eyes "hit" a certain roi
if ~isempty(p.Results.roi)
    factors = horzcat(factors,strcat(p.Results.roi,'_hit'));
end


fixdata = ds(:,horzcat(factors,p.Results.Xvar,p.Results.Yvar));


%if we only want to visualize the hits and not the misses
if p.Results.hitonly && ~isempty(p.Results.roi)
    fixdata = fixdata(fixdata.(strcat(p.Results.roi{1},'_hit'))==1,:);
end


%don't plot when any factor is NaN
if p.Results.skipnan
    
    fixdata = fixdata(sum(isnan(fixdata{:,factors}),2)==0,:);
end


%screen dimensions
screenx = p.Results.screendims(1);
screeny = p.Results.screendims(2);
xcenter = fix(screenx/2);
ycenter = fix(screeny/2);


[fixdata, levels] = makeSubsets(fixdata,factors);


if isempty(p.Results.overlay)
    roi_overlay = zeros(screeny,screenx);
else
    roi_overlay = p.Results.overlay;
end


if length(fixdata)>40
    warning('This would produce over 40 plots. Only plotting the first 40.');
    numplots =40;
else
    numplots = length(fixdata);
end


for s = 1:numplots
    
    %                 temp = fixdata(fixdata.(obj.subject_var)==obj.subs{s},:);
    
    temp = fixdata{s};
    
    figure('MenuBar','none','NumberTitle','off','color','black');
    
    
    %eliminates empty space from each window
%     set(0,'DefaultAxesLooseInset',[0,0,0,0])
    
    %                 if ~isempty(p.Results.overlay)
    imshow(roi_overlay.*0.5);
    hold on;
    %                 end
    
    
    scatter(temp.(p.Results.Xvar),temp.(p.Results.Yvar),p.Results.size,p.Results.color);
    
    if p.Results.zoom
        set(gca,'XLim',[1,screenx]);
        set(gca,'YLim',[1,screeny]);
    end
    
    

    
    %draw a rectangle around where the screen goes, and
    %crosshairs at the center
    if p.Results.screen
        rectangle('Position',[1,1,screenx,screeny],'EdgeColor','w');
        %                     rectangle('Position',[xcenter-10,ycenter-10,20,20],'FaceColor','w','Curvature',[1,1]);
        
        %# vertical line
        hx = graph2d.constantline(xcenter, 'LineStyle',':', 'Color',[.7 .7 .7]);
        changedependvar(hx,'x');
        %# horizontal line
        hy = graph2d.constantline(ycenter,'LineStyle',':', 'Color',[.7 .7 .7]);
        changedependvar(hy,'y');
        
        
    end
    
    %make 0,0 at upper-left so it matches psychtoolbox
    set(gca, 'YDir', 'reverse');
    iptsetpref('ImshowBorder','tight')
    title(fix_title_string(levels.label{s}),'color','white');
end

placeFigures('square',true,'intergapH',10,'intergapV',10,'lowergap',100,'uppergap',25,'toolsize',15);

end