function [  ] = paint( ts, angles, title0, name )
%UNTITLED Summary of this function goes here
set(0,'defaultfigurecolor','w')
plot(ts, angles(:,1), '.r');
hold on;
plot(ts, angles(:,2), '.g');
plot(ts, angles(:,3), '.b');
hold off;
legend('ºá¹ö½Ç','¸©Ñö½Ç', 'º½Ïò½Ç');
title(title0);
a = size(ts);
set(gca, 'XLim', [ts(1) ts(a(2))]);

set(gcf,'position',[200,200,800,400])
F=getframe(gcf);
filename = sprintf('%s.png', name);
imwrite(F.cdata,['./',filename])

end

