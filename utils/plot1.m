hold on
ax = gca;
a = min([ax.XLim(1) ax.YLim(1)]);
b = max([ax.XLim(2) ax.YLim(2)]);
plot([a b], [a b]);