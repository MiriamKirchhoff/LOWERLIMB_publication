function [linear_models] = calculate_main_analysis(all_trials_rm, filter_types)
% function [linear_models] = calculate_main_analysis(all_trials_rm)
%
%
% INPUTS
    % all_trials_rm: Struct with fields:
        % phase:        phase of trial   struct with fields 
            % filter_types ->  double (1 x N_trials)
        % mep:          mep of trial     cell (N_channels x 1)
            % mep types ->                  double (1 x N_trials)
        % bandpower:    bp of trial
            % filter_types ->  double (1 x N_trials)
    % filter_types: Cell with names of filter types
%
% OUTPUTS
    % linear_models: linear model with fields filter_types -> mep types ->
        % Coefficients: Coefficients of fitted model
        % pval:         pvalues, uncorrected
        % R2ordinary   
        % R2adjusted
        % fit:          fitted model
        % pval_corrected:   pvalues corrected over all tests
% 
% version   28.11.2023
% author    Miriam Kirchhoff
% project   C2B

all_pvals = [];
mdl_names = {'phase_cos','phase_sin', 'delta','theta','mu', ...
    'beta', 'gamma_low', 'gamma_high', 'mep'};
idx_mep = 1;

for idx_filter = 1:length(filter_types)

    fprintf(['\n \n Model for ' filter_types{idx_filter} '\n'])

    % select data
    mep = double(all_trials_rm.mep);
    phase = double(all_trials_rm.phase.(filter_types{idx_filter}));
    bp = double(all_trials_rm.bandpower.(filter_types{idx_filter}));

    % fit linear model
    mdl_linear.fit = fitlm(bp, mep, 'VarNames', mdl_names(3:end));

    T = array2table([cos(phase) sin(phase), bp, mep]);
    % Default heading for the columns will be A1, A2 and so on.
    % You can assign the specific headings to your table in the following manner
    T.Properties.VariableNames(:) = mdl_names;

    % fit linear mixed model
    fit_phase = fitlme(T, 'mep ~ phase_cos + phase_sin + delta + theta + mu + beta + gamma_low + gamma_high');
    fit_phase_mu = fitlme(T, 'mep ~ phase_cos + phase_sin + delta + theta + mu + beta + gamma_low + gamma_high + (phase_cos + phase_sin):mu');
    fit_nophase = fitlme(T, 'mep ~ delta + theta + mu + beta + gamma_low + gamma_high');

    % likelihood ratio test
    disp('phase main effect')
    phase_results.(filter_types{idx_filter}) = compare(fit_nophase, fit_phase);
    disp(phase_results.(filter_types{idx_filter}))

    disp('phase mu interaction effect')
    phase_mu_results.(filter_types{idx_filter}) = compare(fit_phase, fit_phase_mu);
    disp(phase_mu_results.(filter_types{idx_filter}))

    mdl_linear.pval = mdl_linear.fit.coefTest;
    mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
    mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;

    temp = dataset2cell(phase_results.(filter_types{idx_filter}));
    temp2 = dataset2cell(phase_mu_results.(filter_types{idx_filter}));
    all_pvals = [all_pvals; mdl_linear.pval; mdl_linear.fit.Coefficients.pValue; temp{3,8}; temp2{3,8}];

    disp(mdl_linear.fit.Coefficients)
    disp(['DFE:' num2str(mdl_linear.fit.DFE)])
    disp(mdl_linear.fit.ModelFitVsNullModel)
    disp(mdl_linear.pval)
    disp(mdl_linear.fit.Rsquared)
    disp(mdl_linear)

    % save model
    linear_models.(filter_types{idx_filter}) = mdl_linear;


end % for idx_filter


%% Multiple testing correction using benjamini-hochberg

reshape_all_pvals = reshape(all_pvals, [], length(filter_types));
multicmp_pvalues = reshape(multicmp(all_pvals,'fdr',0.05), [], length(filter_types));


for idx_filter = 1:length(filter_types)
    linear_models.(filter_types{idx_filter}).pvalue_corrected = ...
        multicmp_pvalues(2:end, idx_filter + length(idx_filter)*(idx_mep-1));
    linear_models.(filter_types{idx_filter}).main_effect_corrected = ...
        multicmp_pvalues(1, idx_filter + length(idx_filter)*(idx_mep-1));

end % for idx_filter


%% Write results into excel spreadsheet

filename = ['LOWERLIMB_publication\results\main_results.xlsx'];

for idx_filter = 1:length(filter_types)
    current_results = linear_models.(filter_types{idx_filter});
    % create table
    results = current_results.fit.Coefficients;
    % add corrected pvals
    results.pValue_corrected = current_results.pvalue_corrected(1:size(results, 1));
    % add columns of NaN for main effects
    results.R2_ordinary(:) = NaN;
    results.R2_adjusted(:) = NaN;
    % add main results
    main_results = table(NaN, NaN, NaN, current_results.pval, ...
        current_results.main_effect_corrected, ...
        current_results.R2ordinary, current_results.R2adjusted, ...
        'VariableNames', results.Properties.VariableNames, ...
        'RowNames', {'Main effect'});
    % add phase results
    ph_results = table(NaN, NaN, NaN, reshape_all_pvals(end-1, idx_filter + length(idx_filter)*(idx_mep-1)), ...
        multicmp_pvalues(end-1, idx_filter + length(idx_filter)*(idx_mep-1)), ...
        NaN, NaN, ...
        'VariableNames', results.Properties.VariableNames, ...
        'RowNames', {'phase'});
    ph_mu_results = table(NaN, NaN, NaN, reshape_all_pvals(end, idx_filter + length(idx_filter)*(idx_mep-1)), ...
        multicmp_pvalues(end, idx_filter + length(idx_filter)*(idx_mep-1)), ...
        NaN, NaN, ...
        'VariableNames', results.Properties.VariableNames, ...
        'RowNames', {'phase:mu'});

    results = [main_results; results; ph_results; ph_mu_results];

    % save as excel sheet
    writetable(results, filename, 'Sheet', [ ...
        filter_types{idx_filter}], 'WriteRowNames',true)

end % for idx_filter

clear current_results


end % eof