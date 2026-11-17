% Open this file in MATLAB and click Run.
run_current();

function run_current()
code = fileparts(mfilename('fullpath'));
model_code = fullfile(code,'..','..','baseline','code');
saved_folder = pwd;
saved_path = path;
cleanup = onCleanup(@()restore_location(saved_folder,saved_path)); %#ok<NASGU>
cd(model_code);
addpath(model_code,'-begin');
clear params physics solver run_model run_steps plot_results plot_steps plot_model plot_reviewed load_results
p=params();
run_model(p,fullfile(code,'..','results','spring_radius'), 'sensitive', fullfile(model_code,'..','results'), fullfile(code,'..','figure'));
end

function restore_location(folder,saved_path)
cd(folder);
path(saved_path);
end
