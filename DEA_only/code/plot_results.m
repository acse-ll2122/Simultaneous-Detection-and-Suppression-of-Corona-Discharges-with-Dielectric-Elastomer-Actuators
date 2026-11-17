function plot_results(out,figure_out)
% Draw the current pulse, base response and reference-radius differences.
if nargin<1
    out = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','pulse_response');
end
if nargin<2,figure_out=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure');end
if ~isfolder(figure_out),mkdir(figure_out);end
[base,study] = params();
summary = readtable(fullfile(out,'pulse_summary.csv'));
supports = {'spring','gas'};
support_titles = {'(a) Spring','(b) Gas'};
coefficients = {'c_1','c_1^{\prime}'};
reference_mm = 1e3*[base.spring.a_ref_m,base.gas.a_ref_m];
colors = [0.12 0.40 0.68;0.80 0.32 0.13;0.16 0.53 0.39];
fig = figure('Visible','on','Color','w','Position',[80 80 1250 990], ...
    'DefaultAxesTitleFontWeight','normal');
% Keep the export canvas independent of visible-window resizing.
set(fig,'PaperUnits','inches','PaperPositionMode','manual', ...
    'PaperPosition',[0 0 fig.Position(3:4)]/get(groot,'ScreenPixelsPerInch'));
for i = 1:numel(supports)
    [curves,cases] = load_pulse_curves(out,summary,supports{i});
    radii_mm = cases.a_ref_mm';
    window_s = [0,cases.simulation_end_ms(1)*1e-3];
    edges_s = [study.on_s,study.pulse_off_s];
    drive = curves{1};
    ax = subplot(5,2,i,'Parent',fig);
    stairs(ax,1e3*drive.time_s,drive.voltage_V/1e3,'k','LineWidth',1.55);
    xlim(ax,1e3*window_s);
    voltage_min_kV = min(drive.voltage_V)/1e3;
    ylim(ax,[1.12*voltage_min_kV,-0.12*voltage_min_kV]);
    grid(ax,'on');
    ylabel(ax,'Voltage (kV)');
    title(ax,sprintf('Applied voltage (0 to %g ms view)',1e3*window_s(2)));

    [~,baseline] = min(abs(radii_mm-reference_mm(i)));
    case_colors = colors(mod(0:numel(curves)-1,size(colors,1))+1,:);
    case_colors(baseline,:) = [0 0 0];
    labels = compose('%g',radii_mm);
    labels{baseline} = sprintf('%g (base case)',radii_mm(baseline));
    ax = subplot(5,2,[i+2,i+4],'Parent',fig);
    handles = plot_window(ax,curves(baseline),case_colors(baseline,:),window_s,edges_s,1);
    xlabel(ax,'');
    title(ax,sprintf('%s   |   %s = %g N s/m^3', ...
        support_titles{i},coefficients{i},cases.damping_N_s_m3(1)));
    key = legend(ax,handles,labels(baseline),'Location','northeast','FontSize',10);
    title(key,'a_{ref} (mm)');

    baseline_um = curves{baseline}.displacement_um;
    for j = 1:numel(curves)
        curves{j}.displacement_um = curves{j}.displacement_um-baseline_um;
    end
    ax = subplot(5,2,[i+6,i+8],'Parent',fig);
    handles = plot_window(ax,curves,case_colors,window_s,edges_s,baseline);
    key = legend(ax,handles,labels,'Location','southeast','FontSize',10);
    title(key,'a_{ref} (mm)');
    limits = ylim(ax);
    padding = 0.04*diff(limits);
    ylim(ax,[min(limits(1),-padding),max(limits(2),padding)]);
    title(ax,'Effect of reference radius over the full window','FontWeight','normal');
    ylabel(ax,sprintf('Difference vs. %g mm (\\mum)',radii_mm(baseline)));
end
sgtitle(fig,sprintf('Figure 3. DEA response to a %g kV pulse from %g to %g ms', ...
    base.pulse_V/1e3,1e3*study.on_s,1e3*study.pulse_off_s));
stem = 'figure_3_pulse_full_difference';
axes_handles = findall(fig,'Type','axes');
for k = 1:numel(axes_handles)
    axes_handles(k).Toolbar.Visible = 'off';
end
drawnow;
% Fixed normalized panel positions preserve the reference Figure 3 layout.
% Order: gas difference, gas response, gas drive, spring difference,
% spring response, spring drive (the reverse axes creation order).
panel_positions = [0.57034090909090907 0.10791336211410715 0.33465909090909096 0.29135776388974455; ...
    0.57034090909090907 0.44670145966032176 0.33465909090909096 0.2913577638897446; ...
    0.57034090909090907 0.78548955720653635 0.33465909090909096 0.12196371511663728; ...
    0.13 0.10791336211410715 0.33465909090909091 0.29135776388974455; ...
    0.13 0.44670145966032176 0.33465909090909091 0.2913577638897446; ...
    0.13 0.78548955720653635 0.33465909090909091 0.12196371511663728];
for k=1:numel(axes_handles)
    set(axes_handles(k),'Units','normalized','Position',panel_positions(k,:));
end
drawnow;
print(fig,fullfile(figure_out,[stem '.png']),'-dpng','-r180');
set(fig,'Visible','on');
savefig(fig,fullfile(figure_out,[stem '.fig']));
% Keep the exported figure open for interactive inspection.
end

function lines = plot_window(ax,curves,colors,window_s,switch_s,emphasis)
% Plot the saved displacement samples and mark the voltage changes.
if nargin<6
    emphasis = [];
end
hold(ax,'on');
lines = gobjects(1,numel(curves));
order = [setdiff(1:numel(curves),emphasis,'stable'),emphasis];
for j = order
    T = curves{j};
    show = T.time_s>=window_s(1)-1e-12 & T.time_s<=window_s(2)+1e-12;
    lines(j) = plot(ax,1e3*T.time_s(show),T.displacement_um(show), ...
        'Color',colors(j,:),'LineWidth',1.15+0.9*ismember(j,emphasis));
end
for t = switch_s
    xline(ax,1e3*t,':','Color',[0.4 0.4 0.4],'HandleVisibility','off');
end
xlim(ax,1e3*window_s);
grid(ax,'on');
box(ax,'on');
set(ax,'FontSize',11,'Layer','top');
xlabel(ax,'Time (ms)');
ylabel(ax,'Radial displacement (\mum)');
end

function [curves,cases] = load_pulse_curves(out,summary,support)
% Load the saved trajectories in the same order as the case summary.
cases = summary(strcmp(summary.support,support),:);
curves = cell(1,height(cases));
for j = 1:height(cases)
    label = strrep(sprintf('%g',cases.a_ref_mm(j)),'.','p');
    file = sprintf('pulse_%s_aref_%smm.csv',support,label);
    curves{j} = readtable(fullfile(out,file));
end
end
