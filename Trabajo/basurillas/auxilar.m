h0 = @(x,y) exp(-((x-0.5).^2+(y-0.5).^2)/0.01);
xi = p(1,:);
yi = p(2,:);
elem = t(1:3,:)';

trisurf(elem,xi,yi,h0(xi,yi))