function hits_all = calcHits(x_coords,y_coords,mask,varargin)
%given a set of x and y coordinates, and a binary mask, determines whether
%each of the coordinates intersects with the mask. Returns a binary vector
%indicating whether or not it "hit" the mask. Can optionally supply a
%screen resolution ('xres',1024,'yres',768), but it assumes it will be the same as your mask. 

p = inputParser;
p.addParameter('xres', size(mask,2),@isnumeric);
p.addParameter('yres', size(mask,1),@isnumeric);
parse(p,varargin{:});

xres = p.Results.xres;
yres = p.Results.yres;

roi_mask = mask;

coords = horzcat(x_coords,y_coords);

%bad samples recoded as zero
coords(coords(:,1)>xres | coords(:,1)<=0 | isnan(coords(:,1)),1) = xres;
coords(coords(:,2)>yres | coords(:,2)<=0 | isnan(coords(:,2)),2) = yres;

coords = ceil(coords);

%find all unique fixations
[ucoords,~,idx_a] = unique(coords,'rows','stable');

scr = zeros(yres,xres);

%find the indicies of the fixations
idx = sub2ind([yres,xres],ucoords(:,2),ucoords(:,1));
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

end
