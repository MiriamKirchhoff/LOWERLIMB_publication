function fig = plot_bp_correlation(data, colors, filter_types)
% function fig = plot_bp_correlation(data, colors)
%
% Plots the correlation of the bandpowers with each other
%
% INPUTS
    % all_trials_rm: Struct with fields:
        % bandpower:        nested struct containing bandpower data
            % C3_C4_Hjorth: bp values for C3_C4 for frequency bands   
            %               double (N_trials x N_bands)
            % C1_C2_Hjorth: bp values for C1_C2 for frequency bands   
            %               double (N_trials x N_bands)
            % names:        frequency band names    cell (1 x 6) 
    % filter_types: Cell with names of filter types
    % colors:  Contains RGB values with fields:
        % darkblue                                  double (1 x 3)
    % filter_types: Cell array of strings containing fieldnames of filters
    %                                               cell (1 x 2)
    % filter_names: Cell array of strings containing titles of filters
    %                                               cell (2 x 1)
%
% OUTPUTS
    % fig:      Figure
% 
% version   28.11.2023
% author    Miriam Kirchhoff
% project   C2B

band_names = {'\delta', '\theta', '\mu', '\beta', '\gamma low', '\gamma high'};
colormap_custom = zeros(100, 3);
for i = 1:3
    colormap_custom(:, i) = linspace(1, colors.darkblue(i), 100);
end

fig = figure('units','normalized','outerposition',[0.6 0.65 0.3 0.35]); 

t = tiledlayout(1, length(filter_types), 'TileSpacing','compact');

title(t, 'Self-correlation of bandpower', 'FontName', 'Times New Roman')

for idx_filter = 1:length(filter_types)
    nexttile
    data_current = round(corr(data.bandpower.(filter_types{idx_filter})), 2);
    h = heatmap(data_current, 'Colormap', colormap_custom, 'FontName', 'TimesNewRoman', 'XDisplayLabels', band_names, 'YDisplayLabels', band_names, 'Colorbarvisible', 'off');

    switch filter_types{idx_filter}
        case 'C3_C4_Hjorth'
            title_text = 'C3/C4 Hjorth';
        case 'C1_C2_Hjorth'
            title_text = 'C1/C2 Hjorth';

    end % switch filter_types{idx_filter}
    title(title_text);
end
fontsize(gcf,scale=1.4)

end % eof