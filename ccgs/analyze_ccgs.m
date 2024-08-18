function ccg_output = analyze_ccgs(data)
% INPUT:
%   data: trials x time x cells

%flag.start_time_index = find(tc>500,1,'first');
%flag.end_time_index = find(tc<1400,1,'last'); %1400
flag.jit_window = 25;
flag.min_lag = -100;
flag.max_lag = 100;


% convert spike train to 'has at least 1 spike in bin'
data.spiketrain = logical(data.spiketrain); 

% get times of interest
%data_subset = data.spiketrain(:, flag.start_time_index:flag.end_time_index, :);
ccg_output = struct;

% jitter spiketrains
for j = 1:size(data_subset, 3)
    data_subset_jitter{j} = jitter(data_subset(:,:,j), flag.jit_window);
    data_subset_real{j} = data_subset(:,:,j);
end

% change to for if do not want to use multiple cpus
Ncell_included=size(data_subset,3);
ccgs=cell(1,Ncell_included);
for j = 1:Ncell_included
    disp("processing neuron: " + j);
    ccgs{j}=struct();
    cnt = 0;
    for k = j+1:Ncell_included
        % pre & post comps
        % first neuron is pre
        cnt = cnt + 1;
        [ccgs{j}.ccg_norm(cnt,:), ccgs{j}.ccg_unnorm(cnt,:), ccgs{j}.gm_FR(cnt,:)] = xcorr_gm(data_subset_real{j}, data_subset_real{k}, flag.max_lag, flag.min_lag);
        [ccgs{j}.ccg_norm_jitter(cnt,:), ccgs{j}.ccg_unnorm_jitter(cnt,:), ccgs{j}.gm_FR_jitter(cnt,:)]  = xcorr_gm(data_subset_jitter{j}, data_subset_jitter{k}, flag.max_lag, flag.min_lag);
        ccgs{j}.neuron_id_pairs(cnt,:) = [j,k];
    end
end
% aggregate ccgs into ccg output struct
fields = fieldnames(ccgs{1});
for j = 1:length(fields)
    ccg_output.(fields{j}) = ccgs{1}.(fields{j});
end

for i = 2:size(data_subset,3)-1
    for j = 1:length(fields)
        ccg_output.(fields{j}) = [ccg_output.(fields{j}); ccgs{i}.(fields{j})];
    end
end

% ccg control is base - jitter
ccg_output.ccg_control = ccg_output.ccg_norm-ccg_output.ccg_norm_jitter;            %normalized,   jitter-corrected CCG (you probably want this)
ccg_output.ccg_control_unnorm = ccg_output.ccg_unnorm-ccg_output.ccg_unnorm_jitter; %unnormalized, jitter-corrected CCG

end