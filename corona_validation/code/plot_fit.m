function plot_fit(U,I,fits,out)
fig=figure('Visible','off','Color','w','Position',[50 50 1050 650]);
set(fig,'PaperUnits','inches','PaperPositionMode','manual', ...
 'PaperPosition',[0 0 fig.Position(3:4)]/get(groot,'ScreenPixelsPerInch'));
ax=axes('Parent',fig);hold(ax,'on');
h(1)=plot(ax,U,I,'ko','MarkerFaceColor','w','LineWidth',1.5);
grid_u=linspace(min(U),max(U),400)';chosen=[1 2 4];colors=lines(3);
for k=1:numel(chosen)
 f=fits(chosen(k));h(k+1)=plot(ax,grid_u,f.K_uA_cm_kVn*max(grid_u-f.U0_kV,0).^f.n, ...
  'Color',colors(k,:),'LineWidth',1.7);
end
xlabel(ax,'Voltage U (kV)');ylabel(ax,'Line current I_l (\muA/cm)');grid(ax,'on');box(ax,'on');
title(ax,'Negative corona: digitized paper data and fitted laws');
legend(ax,h,{'Fig. 3 digitized points','n=2, U_0=21.6 kV, current LS', ...
 'n=2, sqrt-current LS','Free n, fixed onset, log LS'},'Location','northwest');
set(ax,'FontSize',12);drawnow;
print(fig,fullfile(out,'fit_comparison.png'),'-dpng','-r180');
savefig(fig,fullfile(out,'fit_comparison.fig'));close(fig);
end
