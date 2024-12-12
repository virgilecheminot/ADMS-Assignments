function dis_stru(posiz,l,gamma,xy,pr,idb,ndof)
% Plotting the undeformed structure

xmax = max(xy(:,1));
xmin = min(xy(:,1));
ymax = max(xy(:,2));
ymin = min(xy(:,2));

dx = (xmax - xmin)/200;
dy = (ymax - ymin)/200;
d = sqrt(dx^2 + dy^2);

npr     = max(pr);
red     = cos(((0:npr/(npr-1):npr)-npr/2)*pi/2/npr).^2;
green   = cos((0:npr/(npr-1):npr)*pi/2/npr).^2;
blue    = cos(((0:npr/(npr-1):npr)-npr)*pi/2/npr).^2;
colori = [red' green' blue'];

figure();
hold on;
% Step 2: elements
for ii=1:length(posiz)
    xin=posiz(ii,1);
    yin=posiz(ii,2);
    xfi=posiz(ii,1)+l(ii)*cos(gamma(ii));
    yfi=posiz(ii,2)+l(ii)*sin(gamma(ii));
    colore = colori(pr(ii),:);
    plot([xin xfi],[yin yfi],'linewidth',5,'color',colore);
end
grid on; box on;

% Step 1: nodal positions
% plot(xy(:,1),xy(:,2),'r.','markersize',20);

triangolo_h = [ 0 0; -sqrt(3)/2 .5; -sqrt(3)/2 -.5; 0 0]*d*2;
triangolo_v = [ 0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2;
triangolo_r = [0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2 * [sqrt(2)/2 -sqrt(2/2); -sqrt(2)/2 -sqrt(2)/2];


hold on
for ii = 1:size(xy,1)
    rectangle('Position',[xy(ii,1)-d/2 xy(ii,2)-d/2 d d],'curvature',1,'facecolor',[1 0 0],'linewidth',0.1);
    text(xy(ii,1) + d, xy(ii,2) + d,num2str(ii));
    if (idb(ii,1) > ndof)
        fill(xy(ii,1) + triangolo_h(:,1),xy(ii,2) + triangolo_h(:,2),'k');
    end
    if (idb(ii,2) > ndof)
        fill(xy(ii,1) + triangolo_v(:,1),xy(ii,2) + triangolo_v(:,2),'k');
    end
    if (idb(ii,3) > ndof)
        fill(xy(ii,1) + triangolo_r(:,1),xy(ii,2) + triangolo_r(:,2),'k');
    end
end


axis equal
xlim([min(xy(:,1))-10*d max(xy(:,1))+10*d]);
ylim([min(xy(:,2))-10*d max(xy(:,2))+10*d]);


title('Undeformed Structure')