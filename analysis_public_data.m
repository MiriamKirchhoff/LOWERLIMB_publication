%% LOWERLIMB: Data analysis file for publication
%
% This file includes the data analysis and figures that do not require the
% raw data. For privacy reasons, please contact Ulf Ziemann for access to
% the raw data set.
%
% author    Miriam Kirchhoff
% project   C2B

%% Load and restructure data

data = read_and_restructure_data();


%% Get settings

filter_types = fieldnames(data.phase);

% define names
bp_names = {'\delta', '\theta', '\mu', '\beta', '\gamma low', '\gamma high'};
filter_names = {'C3/C4 Hjorth', 'C1/C2 Hjorth'};

% Define colors
colors.blue     = [25 200 255]/255;
colors.yellow  = [153 140 28]/255;
colors.red      = [255 25 86]/255;
colors.darkred  = [143 43 67]/255;
colors.darkblue = [57 86 138]/255;
colors.grey     = [1 1 1]*0.75;


%% Prepare dataset: Exclude outlier trials, log-transform, and z-score

% Exclude outlier trials
fields = fieldnames(data);
for field_idx = 1:length(fields)
    field_name = fields{field_idx};
    temp = data.(field_name);

    % iterate over fields for nested structs, except for outlier field
    if isstruct(data.(field_name)) & ~ matches(field_name, 'outlier')

        data_rm.(field_name).C3_C4_Hjorth = temp.C3_C4_Hjorth(~data.outlier.all,:);
        data_rm.(field_name).C1_C2_Hjorth = temp.C1_C2_Hjorth(~data.outlier.all,:);

    elseif ~ isstruct(data.(field_name))
        data_rm.(field_name) = temp(~data.outlier.all,:);

    end
end

% Z-score participant-wise
participants = unique(data_rm.participant);
for idx_participant = 1:length(participants)
    participant_location = data_rm.participant == participants(idx_participant);

    data_rm.mep(participant_location) = zscore(log(data_rm.mep(participant_location)));
    data_rm.bandpower.C3_C4_Hjorth(participant_location,:) = ...
        zscore(log(data_rm.bandpower.C3_C4_Hjorth(participant_location,:)));
    data_rm.bandpower.C1_C2_Hjorth(participant_location,:) = ...
        zscore(log(data_rm.bandpower.C1_C2_Hjorth(participant_location,:)));

    data_rm.bandpower_all.C3_C4_Hjorth(participant_location,:) = ...
        zscore(log(data_rm.bandpower_all.C3_C4_Hjorth(participant_location,:)));
    data_rm.bandpower_all.C1_C2_Hjorth(participant_location,:) = ...
        zscore(log(data_rm.bandpower_all.C1_C2_Hjorth(participant_location,:)));
end

data_rm.bandpower_all.frequency = data.bandpower_all.frequency;

clear temp field_idx idx_participant participants participant_location ...
    fields field_name data

%% Calculate effect of bandpower on MEPs

bp_stats = [];

for idx_filter = 1:length(filter_types)
    for idx_bp = 1:size(data_rm.bandpower_all.(filter_types{idx_filter}),2)

        % select current data
        mep = data_rm.mep;
        phase = data_rm.phase.(filter_types{idx_filter});
        bp = data_rm.bandpower_all.(filter_types{idx_filter})(:,idx_bp);

        % calculate correlation
        [correlation, corr_p, corr_CI_low, corr_CI_high] = corrcoef(bp', mep');

        bp_stats.(filter_types{idx_filter}).correlation(idx_bp,:) = correlation(2,1);
        bp_stats.(filter_types{idx_filter}).correlation_CI(idx_bp,:) = [corr_CI_low(2,1); corr_CI_high(2,1)];

    end  % for idx_bp
end % for idx_mep

clear mep bp phase idx_filter idx_bp correlation corr_p corr_CI_low corr_CI_high


%% Calculate effect of pre-defined frequency bands

disp('Start calculating main analysis'); tic

% Calculate and save main effects
main_results = calculate_main_analysis(data_rm, filter_types);

fprintf('Completed calculating and saving main analysis. Elapsed time: %.0f seconds \n', toc)


%% Calculate post-hoc analysis: Interaction with phase

disp('Start calculating main analysis'); tic

% Calculate and save main effects
post_hoc_results = calculate_phase_interaction_analysis(data_rm, filter_types);

fprintf('Completed calculating and saving main analysis. Elapsed time: %.0f seconds \n', toc)


%% Plot effect of bandpower on MEPs

disp('Start plotting dependence of bandpower and MEPs'); tic

fig = plot_bp_mep_dependence(data_rm, bp_stats, main_results, colors, filter_types, filter_names, bp_names);
saveas(fig, 'LOWERLIMB_publication\figures\bp_mep_dependence.png')

fprintf('Completed + saved plot. Elapsed time: %.0f seconds \n', toc)
clear fig


%% Scatterplot of each bandpower and MEP and phase and MEP

disp('Start plotting dependence of phase-bandpower and MEPs'); tic

fig = plot_bp_phase_mep(data_rm, colors, filter_types, filter_names);
saveas(fig, 'LOWERLIMB_publication\figures\Correlation_bp_phase.jpg')

fprintf('Completed + saved plot. Elapsed time: %.0f seconds \n', toc)
clear fig


%% Plot correlation of bandpower with each other

disp('Start plotting dependence of phase-bandpower and MEPs'); tic

fig = plot_bp_correlation(data_rm, colors, filter_types);
saveas(fig,'LOWERLIMB_publication\figures\bandpowercorrelation.jpg')

fprintf('Completed + saved plot. Elapsed time: %.0f seconds \n', toc)
clear fig



