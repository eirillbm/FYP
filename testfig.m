figure(1);

box off;

set(gca, 'FontName','Times New Roman','FontSize', 12);

xlabel('$${x_{12}*}$$','Interpreter','Latex');

ylabel('Probability density');

set(legend,'FontName','Times New Roman','FontSize',8);

 legend boxoff;

 h=gcf;

set(gcf,'PaperPosition', [0 0 10 7]);

set(gcf,'PaperSize', [10 7]);

ext = '.pdf';

print(gcf, '-dpdf','numberofwaves');