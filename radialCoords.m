function coords = radialCoords(xcenter,ycenter,numpoints,radius,starting_angle)

coords=zeros(numpoints,2);

for i=1:numpoints
angle=(i-1)*(360/numpoints)+starting_angle;
coords(i,1)=ceil(xcenter+sin((angle*pi)/180)*radius);
coords(i,2)=ceil(ycenter+cos((angle*pi)/180)*radius);

end


