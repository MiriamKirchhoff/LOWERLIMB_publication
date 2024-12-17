function [fig] = plot_bp_mep_dependence(data, bp_stats, linear_models, colors, filter_types, filter_names, bp_names)
% function [fig] = plot_bp_mep_dependence(data, bp_stats, linear_models, colors)
%
% Plot dependence of MEPs on bandpowers using correlation, marking
% significant relationships of frequency bands with MEPs
%
% INPUTS
    % data: Struct with fields:
        % phase:        phase of trial   struct with fields 
            % filter_types -> mep types ->  double (1 x N_trials)
        % mep:          mep of trial        cell (N_channels x 1)
            % mep types ->                  double (1 x N_trials)
        % bandpower:    bp of trial
            % filter_types -> mep types ->  double (1 x N_trials)
    % bp_stats: Struct with fields filter_types -> mep types ->
        % correlation   corr of bp and MEP  
    % linear_models: linear model with fields filter_types -> mep types ->
        % Coefficients: Coefficients of fitted model
        % pval:         pvalues, uncorrected double(1 x N_frequencybands)
        % R2ordinary   
        % R2adjusted
        % fit:          fitted model
        % pval_corrected:   pvalues corrected over all tests
    % freq_band:    Name of the frequency band  str
    % colors:   Struct with fields:
        % red           rgb color vector    double (1x3)
        % darkred       rgb color vector    double (1x3)
    % mep_types:    Cell with names of mep types
    % filter_types: Cell with names of filter types
%
% OUTPUTS
    % fig           figure handle       fig
% 
% version   1.0, 28.11.2023
% author    Miriam Kirchhoff
% project   C2B


%% Init

% Pre-defined bandpower frequencies in [Hz]
bandpowerfreq = ...
    [0.5 3;     % delta
    4 7;        % theta
    8 12.5;     % mu
    13 30;      % beta
    30.5 45;    % gamma low
    55 95];     % gamma high

% Discreticized bandpowers in [Hz]
bandpowers = data.bandpower_all.frequency;


%% Plotting

fig = figure('units','normalized','outerposition',[0.6 0.4 0.4 0.6]); 
set(gcf,'color','w');
t = tiledlayout(3,6, 'TileSpacing','tight');
title(t, 'Correlation of bandpower of different frequency and MEPs', 'FontName', 'Times New Roman')
lim = [-0.07 0.07];


for idx_filter = 1:length(filter_types)

    nexttile([1,3]); hold on; grid on; set(gca, 'FontName', 'Times New Roman')

    % plot mean correlation
    plot(bandpowers, bp_stats.(filter_types{idx_filter}).correlation, '-', 'color', 'k', 'LineWidth',2)

    % plot CI of correlation
    X_area = [bandpowers', flip(bandpowers')];
    ci = bp_stats.(filter_types{idx_filter}).correlation_CI;
    Y = [ci(:,1), flip(ci(:,2))];
    fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)

    % Plot 50Hz noise markers
    rectangle('Position',[45 lim(1) 10 lim(2)-lim(1)], 'FaceColor', [1 1 1 0.9], 'EdgeColor', [1 1 1 0.9])
    rectangle('Position',[95 lim(1) 10 lim(2)-lim(1)], 'FaceColor', [1 1 1 0.9], 'EdgeColor', [1 1 1 0.9])

    mdl_linear = linear_models.(filter_types{idx_filter});

    % find significant relationships
    significance_idx = mdl_linear.pvalue_corrected < 0.05;


    for i = 2:size(bandpowerfreq,1) + 1
        col = [1 1 1]*0.4;
        % re-color for significant
        if significance_idx(i)
            col = colors.red;
        end

        % draw box
        % rectangle('Position',[bandpowerfreq(i-3,1) lim(1) bandpowerfreq(i-3,2)-bandpowerfreq(i-3,1) lim(2)-lim(1)], 'FaceColor', [colors.blue 0.05], 'EdgeColor', [colors.blue 0.05])

        % indicate slope
        stem(mean(bandpowerfreq(i-1,:)), mdl_linear.fit.Coefficients.Estimate(i), 'filled', 'LineWidth',2, 'Color', col, 'Marker','.', 'MarkerSize',20)
        beta = mdl_linear.fit.Coefficients.Estimate(i);
        %scatter(mean(bandpowerfreq(i-1,:)), beta + 0.01*sign(beta), 60, 'wo', 'filled')

        %scatter(mean(bandpowerfreq(i-1,:)), beta + 0.01*sign(beta), 60, 'k*')
    end

    xline(bandpowerfreq(:,1), '-', ...
        {'\delta', '\theta', '\mu', '\beta', '\gamma low', '\gamma high'}, 'LabelOrientation', 'horizontal', 'FontName', 'Times New Roman')
    yline(0)
    xlabel('Bandpower frequency')
    if idx_filter == 1; ylabel('Correlation');
    else; yticklabels([]); end
    xticks(sort([bandpowerfreq(2:end,1); 45; 95; 10; 20; 30; 40; 50; 60; 70; 80; 90]))
    xticklabels(["4", "8", "", "13", "20", "", "30.5", "", "45", "", "55", "", "70", "80", "", "95"])
    xtickangle(0)
    yticks(-0.1:0.02:0.1)
    title_text = '';



    if length(filter_types) > 1
        switch filter_types{idx_filter}
            case 'SSD'
                title_text = [title_text 'SSD'];
            case 'C3_4_Hjorth'
                title_text = [title_text 'C3/C4 Hjorth'];
            case 'C1_2_Hjorth'
                title_text = [title_text 'C1/C2 Hjorth'];

        end % switch filter_types{idx_filter}
    end % if length(filter_types) > 1

    title(title_text);
    xlim([0 95]); ylim(lim)

end % for idx_filter

lgd = legend("Correlation \pm 95% CI", "", "Estimate of significant effects", 'location', 'north', 'numcolumns', 2);
lgd.Layout.Tile = 'north';

for idx_filter = 1:length(filter_types)
    for idx_bandpower = 1:6
        nexttile

        current_bp = data.bandpower.(filter_types{idx_filter})(:,idx_bandpower);
        scatter(current_bp, data.mep, 2, 'k', 'filled', MarkerFaceAlpha= 0.2, MarkerEdgeAlpha= 0.2)
        current_correlation = linear_models.(filter_types{idx_filter}).fit.Coefficients.Estimate(idx_bandpower + 1);
        current_pval = linear_models.(filter_types{idx_filter}).pvalue_corrected(idx_bandpower +1);
        x = [min(current_bp),max(current_bp)];
        y = x*current_correlation;

        hold on

        if current_pval > 0.05
            plot(x,y, LineWidth=2, Color=[1 1 1]*0.8)
        else
            plot(x,y, LineWidth=2, Color=colors.red)
        end

        grid on

        xlim([-3 3])
        ylim([-3 3])
        axis square
        if idx_bandpower == 1
            % give ylabel
            ylabel({['\bf{' filter_names{idx_filter} '}'];'\rm MEP (log, z-score)'})
        else
            yticklabels([])

        end
        if idx_filter == 2
            xlabel([bp_names{idx_bandpower} ' power'])
        else
            xticklabels([])
        end
    end

end

xlabel(t, {'', 'Band power (z-score)'}, 'FontName', 'Times New Roman', 'FontSize', 9)

fontsize(gcf,scale=1.4)
