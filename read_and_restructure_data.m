function data = read_and_restructure_data()
% function data = read_and_restructure_data()
%
%
%
% OUTPUTS
    % data: struct containing the restructured data, with fields:
        % mep:              mep of each trial       double (N_trials x 1)
        % participant:      particiant number       double (N_trials x 1)
        % preinnervation:   Preinnervation values   double (N_trials x 1)
        % outlier:          nested struct containing outlier labels
            % preinnervation                        double (N_trials x 1)
            % EEG_threshold                         double (N_trials x 1)
            % all                                   logical (N_trials x 1)
        % phase:            nested struct containing phase data
            % C3_C4_Hjorth: Phase (C3_C4)           double (N_trials x 1)
            % C1_C2_Hjorth: Phase (C1_C2)           double (N_trials x 1)
        % bandpower:        nested struct containing bandpower data
            % C3_C4_Hjorth: bp values for C3_C4 for frequency bands   
            %               double (N_trials x N_bands)
            % C1_C2_Hjorth: bp values for C1_C2 for frequency bands   
            %               double (N_trials x N_bands)
            % names:        frequency band names    cell (1 x 6) 
        % bandpower_all:    nested struct with continuous bandpower data
            % C3_C4_Hjorth: continuous bp values, C3_C4 across frequencies   
            %               double (N_trials x N_freqs)
            % C1_C2_Hjorth: continuous bp values, C1_C2 across frequencies   
            %               double (N_trials x N_freqs)
            % frequency:    corresponding frequency values 
            %               double (1 x N_freqs)
%
% version   17.12.2024
% author    Miriam Kirchhoff
% project   C2B

%% Load the data from csv

data_loaded = readtable("analysis_data.csv");


%% Re-structure data

data.mep            = data_loaded.mep;
data.participant    = data_loaded.participant;
data.preinnervation = data_loaded.preinnervation;

data.outlier.preinnervation = data_loaded.outlier_preinnervation;
data.outlier.EEG_threshold = data_loaded.outlier_EEG_threshold;
% data.outlier.no_peak = data_loaded.outlier_no_peak;
data.outlier.all = data.outlier.preinnervation | data.outlier.EEG_threshold;

data.phase.C3_C4_Hjorth = data_loaded.phase_C3_4_Hjorth;
data.phase.C1_C2_Hjorth = data_loaded{:,"phase_C1_2_Hjorth"};


%% Bandpower

% bandpower bins
idx_start = find(matches(data_loaded.Properties.VariableNames, 'C3_4_Hjorth_bandpower_delta'));
idx_end = find(matches(data_loaded.Properties.VariableNames, 'C3_4_Hjorth_bandpower_gamma_high'));
data.bandpower.C3_C4_Hjorth = data_loaded{:, idx_start:idx_end};
idx_start = find(matches(data_loaded.Properties.VariableNames, 'C1_2_Hjorth_bandpower_delta'));
idx_end = find(matches(data_loaded.Properties.VariableNames, 'C1_2_Hjorth_bandpower_gamma_high'));
data.bandpower.C1_C2_Hjorth = data_loaded{:, idx_start:idx_end};
data.bandpower.names = {'delta', 'theta', 'mu', 'beta', 'gamma_low', 'gamma_high'};

% bandpower continuous
idx = contains(data_loaded.Properties.VariableNames, 'C3_4_Hjorth_bandpower_cont');
data.bandpower_all.C3_C4_Hjorth = data_loaded{:, idx};
idx = contains(data_loaded.Properties.VariableNames, 'C1_2_Hjorth_bandpower_cont');
data.bandpower_all.C1_C2_Hjorth = data_loaded{:, idx};

% get frequencies of continuous bandpower measures from header
idx = find(idx);
bp_value = ones(size(idx));
for i = 1:length(idx)
    temp = split(data_loaded.Properties.VariableNames(idx(i)), '_');
    if length(temp) == 6
        bp_value(i) = str2double(temp{6});
    else % if length(temp) == 6
        bp_value(i) = str2double(temp{6}) + str2double(['0.' temp{7}]);
    end % if length(temp) == 6
end % for i = 1:length(idx)
data.bandpower_all.frequency = bp_value;

end % eof
