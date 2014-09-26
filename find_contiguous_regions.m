function b=find_contiguous_regions(V,varargin)
%%

if ~isempty(varargin)
    val = varargin{1};
    
    V = V ~= val;
else

V=~isnan(V);
end


D = diff([0,V,0]);
    b.regstart = find(D == 1);
    b.regend = find(D == -1) - 1;
    b.regwidth=b.regend-b.regstart+1;


end