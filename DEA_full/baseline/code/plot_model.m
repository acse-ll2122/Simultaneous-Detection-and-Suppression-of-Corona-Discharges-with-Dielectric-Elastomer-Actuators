function plot_model(out,figure_out)
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,out=fullfile(root,'results');end
if nargin<2,figure_out=fullfile(root,'figure');end
plot_reviewed(out,figure_out,'baseline');
end
