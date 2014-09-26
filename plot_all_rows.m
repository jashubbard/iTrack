function h=plot_all_rows(data,varargin)





if nargin>1
    h=varargin{1};
    subplot(h);
else
    h=figure;
end

% cmap=colormap(hsv(size(data,1)));
cmap=cbrewer('qual','Set1',size(data,1));

for i=1:size(data,1)
    y=data(i,:);
    x=1:size(y,2);
    
    plot(x,y,'Linewidth',2,'color',cmap(i,:)); hold all;

end


ylim([min(data(:))-min(data(:))*.20, max(data(:))*1.2]);


hold off;

return