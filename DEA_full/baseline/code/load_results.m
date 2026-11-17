function S = load_results(out)
% Read saved trajectories and resolve a shared baseline within this package.
S = load(fullfile(out,'run.mat'));
if isfield(S,'shared_baseline') && ~isempty(S.shared_baseline)
    reference = load(fullfile(out,S.shared_baseline));
    assert(isequaln(reference.base,S.base),'The shared baseline parameters changed; rerun the sweep.');
    index = find(strcmp(reference.names,S.names{S.baseline_spring}),1);
    assert(~isempty(index),'Shared baseline trajectory not found.');
    S.cases{S.baseline_spring} = reference.cases{index};
end
end
