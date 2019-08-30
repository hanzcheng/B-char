
x=2:2:50;

y=max(errG(:,2:2:50));

plot(x,y,'b-o');

xl = get(gca,'XLabel');
xlFontSize = get(xl,'FontSize');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 24)
set(xl, 'FontSize', 24);


xlabel('N');

yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 24)
set(yl, 'FontSize', 24);

ylabel('Error in mass conservation');