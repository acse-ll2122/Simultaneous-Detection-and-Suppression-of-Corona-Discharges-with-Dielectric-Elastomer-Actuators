function plot_steps(out,figure_out)
if nargin<1,out=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','step_damping');end
if nargin<2,figure_out=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure');end
if ~isfolder(figure_out),mkdir(figure_out);end
S=load(fullfile(out,'step.mat'));colors=lines(5);colors(3,:)=[0 0 0];
fig=figure('Visible','off','Color','w','Position',[50 50 1320 1128], ...
 'DefaultAxesTitleFontWeight','normal');
ax=axes('Parent',fig,'Units','normalized','Position',[.09 .78 .86 .12]);
plot(ax,[0 10 10 50],[0 0 -1.5 -1.5],'k-','LineWidth',1.8);
xlim(ax,[0 50]);ylim(ax,[-1.65 .15]);grid(ax,'on');box(ax,'on');
title(ax,'Applied voltage');xlabel(ax,'Time (ms)');ylabel(ax,'Voltage (kV)');set(ax,'FontSize',12);
labels=compose('%g',S.dampings);labels{3}='1000 (base)';
styles={'-','--','-','-.',':'};
names={'Spring','Gas'};
for i=1:2
 for j=1:2
  ax=axes('Parent',fig,'Units','normalized','Position',[.09+.49*(i-1),.09+.36*(2-j),.37,.27]);hold(ax,'on');handles=gobjects(1,5);
  for k=1:5
   T=S.curves{i,k}; handles(k)=plot(ax,1e3*T.time_s,T.displacement_um,'Color',colors(k,:),'LineStyle',styles{k},'LineWidth',1.45);
  end
  if j==1,xlim(ax,[0 50]);title(ax,names{i});else,xlim(ax,[10 13]);title(ax,[names{i} ' onset detail']);end
  xlabel(ax,'Time (ms)');ylabel(ax,'Radial displacement (\mum)');grid(ax,'on');box(ax,'on');set(ax,'FontSize',12);
  if j==1,key=legend(ax,handles,labels,'Location','southeast','NumColumns',2,'FontSize',10);title(key,'Damping (N s/m^3)');end
 end
end
annotation(fig,'textbox',[0 .965 1 .03],'String','Figure 2  DEA step response and damping', ...
 'EdgeColor','none','HorizontalAlignment','center','FontSize',16);
drawnow;
if exist('exportgraphics','file')
 if verLessThan('matlab','25.1'),exportgraphics(fig,fullfile(figure_out,'figure_2_step_damping.png'),'Resolution',180);
 else,exportgraphics(fig,fullfile(figure_out,'figure_2_step_damping.png'),'Resolution',180,'Padding',20);end
else,print(fig,fullfile(figure_out,'figure_2_step_damping.png'),'-dpng','-r180');end
savefig(fig,fullfile(figure_out,'figure_2_step_damping.fig'));close(fig);
end
