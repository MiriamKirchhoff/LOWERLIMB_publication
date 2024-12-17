function [linear_models] = calculate_phase_interaction_analysis(data, filter_types)
% function [linear_models] = calculate_phase_interaction_analysis(data, filter_types);
%
% Calculates post-hoc analysis for phase interaction with bandpowers
%
% INPUTS
    % data: Struct with fields:
        % phase:        phase of trial   struct with fields 
            % filter_types ->  double (1 x N_trials)
        % mep:          mep of trial  double (1 x N_trials)
        % bandpower:    bp of trial
            % filter_types ->  double (1 x N_trials)
    % filter_types: Cell with names of filter types
%
% OUTPUTS
    % linear_models: linear model with fields filter_types 
        % Coefficients: Coefficients of fitted model
        % pval:         pvalues, uncorrected
        % R2ordinary   
        % R2adjusted
        % fit:          fitted model
        % pval_corrected:   pvalues corrected over all tests
% 
% version   11.12.2023
% author    Miriam Kirchhoff
% project   C2B

all_pvals = [];
mdl_names = {'phase_cos','phase_sin', 'delta','theta','mu', ...
    'beta', 'gamma_low', 'gamma_high', 'mep'};

for idx_filter = 1:length(filter_types)

    fprintf(['\n \n Model for ' filter_types{idx_filter} '\n'])

    % select data
    mep = double(data.mep);
    phase = double(data.phase.(filter_types{idx_filter}));
    bp = double(data.bandpower.(filter_types{idx_filter}));

    % fit linear model
    mdl_linear.fit = fitlm([cos(phase) sin(phase), bp], mep, ...
        'VarNames', mdl_names);

    T = array2table([cos(phase) sin(phase), bp, mep]);
    % Default heading for the columns will be A1, A2 and so on.
    % You can assign the specific headings to your table in the following manner
    T.Properties.VariableNames(:) = mdl_names;

    % fit linear mixed model
    fit_phase = fitlme(T, 'mep ~ phase_cos + phase_sin + delta + theta + mu + beta + gamma_low + gamma_high');
    fit_nophase = fitlme(T, 'mep ~ delta + theta + mu + beta + gamma_low + gamma_high');

    interaction_pvals = [];
    for idx_variable = 1:length(mdl_names) - 3 % iterate over variables except dependent variable and

        variable = mdl_names{idx_variable + 2};
        fit_phase_interaction = fitlme(T, ...
            ['mep ~ phase_cos + phase_sin + delta + theta + mu + beta' ...
            '+ gamma_low + gamma_high + (phase_cos + phase_sin):' variable]);

        temp_results = compare(fit_phase, fit_phase_interaction);
        %disp(temp_results)
        temp = dataset2cell(temp_results);
        temp_pval = temp{3, 8};
        interaction_pvals(idx_variable) = temp_pval;
        if temp_pval < 0.05
            fprintf([variable ' interaction potentially significant: p = %0.3f \n'], temp_pval)
        else
            fprintf([variable ' interaction not significant: p = %0.3f \n'], temp_pval)
        end

    end % for idx_variable

    mdl_linear.pval = mdl_linear.fit.coefTest;
    mdl_linear.R2ordinary = mdl_linear.fit.Rsquared.Ordinary;
    mdl_linear.R2adjusted = mdl_linear.fit.Rsquared.Adjusted;
    all_pvals = [all_pvals; mdl_linear.pval; mdl_linear.fit.Coefficients.pValue; interaction_pvals'];

    % save model
    linear_models.(filter_types{idx_filter}) = mdl_linear;

end % for idx_filter


%% Multiple testing correction using benjamini-hochberg

multicmp_pvalues = reshape(multicmp(all_pvals,'fdr',0.05), [], length(filter_types));

for idx_filter = 1:length(filter_types)
    linear_models.(filter_types{idx_filter}).pvalue_corrected = ...
        multicmp_pvalues(2:end, idx_filter);
    linear_models.(filter_types{idx_filter}).main_effect_corrected = ...
        multicmp_pvalues(1, idx_filter);

end % for idx_filter


%% Write results into excel spreadsheet

filename = ['LOWERLIMB_publication\results\post_hoc_results.xlsx'];

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
    temp = table(NaN, NaN, NaN, current_results.pval, ...
        current_results.main_effect_corrected, ...
        current_results.R2ordinary, current_results.R2adjusted, ...
        'VariableNames', results.Properties.VariableNames, ...
        'RowNames', {'Main effect'});
    results = [temp; results];

    % save
    writetable(results, filename, 'Sheet', [filter_types{idx_filter}], ...
        'WriteRowNames',true)

end % for idx_filter

clear current_results


end % eof