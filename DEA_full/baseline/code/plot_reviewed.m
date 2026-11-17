function plot_reviewed(out,figure_out,mode)
% Replot the selected study using saved trajectories only.
if nargin<2,figure_out=fullfile(fileparts(out),'figure');end
if nargin<3,mode='baseline';end
S=load_results(out);
if ~isfolder(figure_out),mkdir(figure_out);end
if strcmp(mode,'baseline'),plot_supports(S,figure_out);
elseif strcmp(mode,'sensitive'),plot_sweep(S,figure_out);
else,error('DEA:PlotMode','Unknown study mode.');end
end

function plot_supports(S,out)
fig=figure('Visible','off','Color','w','Position',[40 40 1360 1240], ...
 'DefaultAxesTitleFontWeight','normal');
for col=1:2
 if col==1,index=S.baseline_spring;finish=5;label='Spring';else,index=numel(S.cases);finish=1;label='Gas';end
 unit='s';factor=1;
 T=window(S.cases{index},finish);t=T.time_s*factor;
 % Retain the source figure's paired mechanical and voltage upper panels.
 ax=panel(fig,col,1);colororder(ax,[.04 .32 .66;.84 .26 .12]);
 yyaxis(ax,'left');
 v=T.radius_m*1000;yl='Radius r (mm)';key='r';
 hm=plot(ax,t,v,'-','LineWidth',1.8);ylabel(ax,yl);
 yyaxis(ax,'right');v=T.V_DEA_V;vl='V_{DEA} (V)';vk='V_{DEA}';
 hv=plot(ax,t,v,'--','LineWidth',1.8);ylabel(ax,vl);
 title(ax,sprintf('%s  0-%g %s',label,finish*factor,unit));
 legend(ax,[hm hv],{key,vk},'Location','best','FontSize',10);
 format_ax(ax,finish*factor,unit);
 % Field, corona current and signed leakage each have their own subplot.
 fields={'E4_V_m','jc_A_m2','jo_A_m2'};
 labels={{'Outer surface field','E_4 (MV/m)'},'Corona j_c (\muA/m^2)','Leakage j_o (\muA/m^2)'};
 scales=[1e-6 1e6 1e6];
 for q=1:3
  ax=panel(fig,col,q+1);cc=[.04 .32 .66];st='-';
  if q==3,cc=[.84 .26 .12];st='--';end
  plot(ax,t,T.(fields{q})*scales(q),st,'Color',cc,'LineWidth',1.8);
  ylabel(ax,labels{q});format_ax(ax,finish*factor,unit);
  if q==1,threshold_line(ax,S.base.Eth_V_m*1e-6,finish,col==2);end
 end
end
name='support_components';tt='Coupled response: Spring 5 s and Gas 1 s';
sgtitle(fig,tt,'FontSize',16);
save_plot(fig,out,name);
end

function ax=panel(fig,col,row)
ax=axes('Parent',fig,'Units','normalized','Position',[.105+.49*(col-1),.095+.22*(4-row),.315,.15]);
end

function plot_sweep(S,out)
fig=figure('Visible','off','Color','w','Position',[40 40 1450 1090]);
fields={'radius_m','V_DEA_V','E4_V_m','jc_A_m2','jo_A_m2'};
units={'Radius r (mm)','V_{DEA} (V)',{'Outer surface field','E_4 (MV/m)'},'Corona j_c (\muA/m^2)','Leakage j_o (\muA/m^2)'};
factors=[1e3 1 1e-6 1e6 1e6];colors=lines(numel(S.spring_mm));colors(S.baseline_spring,:)=[0 0 0];
finish=5;time_scale=1;time_unit='s';name='spring_components';
caption='Coupled Spring response for different free radii over 5 s';
styles={'--','-.','-',':'};markers={'o','s','d','^'};
positions={[.10 .70 .35 .19],[.59 .70 .35 .19],[.10 .40 .84 .19],[.10 .10 .35 .19],[.59 .10 .35 .19]};
for q=1:5
 ax=axes('Parent',fig,'Units','normalized','Position',positions{q});hold(ax,'on');hh=gobjects(1,numel(S.spring_mm));
 for k=1:numel(S.spring_mm)
  T=window(S.cases{k},finish);mt=[.1 .35 .7 1.2 2 3 4 4.8]+.015*(k-1);
  ids=unique(interp1(T.time_s,(1:height(T))',mt,'nearest'));
  values=T.(fields{q})*factors(q);
  hh(k)=plot(ax,T.time_s*time_scale,values,'Color',colors(k,:), ...
   'LineStyle',styles{k},'LineWidth',1.7,'Marker',markers{k},'MarkerIndices',ids,'MarkerSize',5,'MarkerFaceColor','w');
 end
 if q==1
  ll=compose('%g mm',S.spring_mm);ll{S.baseline_spring}='3 mm (base)';
  location='northeast';
  legend(ax,hh,ll,'Location',location,'NumColumns',2,'FontSize',11);
 end
 ylabel(ax,units{q});format_ax(ax,finish*time_scale,time_unit);
 if q==3,threshold_line(ax,S.base.Eth_V_m*1e-6,finish,false);end
end
sgtitle(fig,caption,'FontSize',16);
save_plot(fig,out,name);
end

function T=window(T,finish)
if T.time_s(end)>finish
 last=array2table(interp1(T.time_s,T{:,:},finish),'VariableNames',T.Properties.VariableNames);
 T=[T(T.time_s<finish,:);last];
end
end
function format_ax(ax,finish,unit)
xlim(ax,[0 finish]);xlabel(ax,['Time (' unit ')']);grid(ax,'on');box(ax,'on');set(ax,'FontSize',12);
end
function threshold_line(ax,threshold,finish,label_below)
% Keep the actual threshold visible, including when it was on the lower border.
hold(ax,'on');limits=ylim(ax);
if label_below
 limits(1)=min(limits(1),threshold-0.02);ylim(ax,limits);
 direction=-1;alignment='top';offset=0.025;
else
 direction=1;alignment='bottom';offset=0.08;
end
yline(ax,threshold,'--','Color',[0 0 0],'LineWidth',2, ...
 'HandleVisibility','off','Tag','corona_threshold');
text(ax,0.03*finish,threshold+direction*offset*diff(limits), ...
 sprintf('E_{th} = %g MV/m',threshold),'Interpreter','tex', ...
 'FontSize',12,'FontWeight','bold','Color',[0 0 0], ...
 'BackgroundColor','w','Margin',1,'HorizontalAlignment','left', ...
 'VerticalAlignment',alignment,'Tag','corona_threshold_label');
end
function save_plot(fig,out,name)
set(fig,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 fig.Position(3:4)]/get(groot,'ScreenPixelsPerInch'));
drawnow;
if exist('exportgraphics','file')
 if verLessThan('matlab','25.1'),exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',180);
 else,exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',180,'Padding',20);end
else,print(fig,fullfile(out,[name '.png']),'-dpng','-r180');end
savefig(fig,fullfile(out,[name '.fig']));close(fig);
end
