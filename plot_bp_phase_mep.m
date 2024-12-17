function fig = plot_bp_phase_mep(data, colors, filter_types, filter_names)
% function fig = plot_bp_phase_mep(data, colors, filter_types, filter_names)
% 
% Plots phase-bandpower interaction
%
% INPUTS
    % data: struct containing the restructured data, with fields:
        % mep:              mep of each trial       double (N_trials x 1)
        % phase:            nested struct containing phase data
            % C3_C4_Hjorth: Phase (C3_C4)           double (N_trials x 1)
            % C1_C2_Hjorth: Phase (C1_C2)           double (N_trials x 1)
        % bandpower:        nested struct containing bandpower data
            % C3_C4_Hjorth: bp values for C3_C4 for frequency bands   
            %               double (N_trials x N_bands)
            % C1_C2_Hjorth: bp values for C1_C2 for frequency bands   
            %               double (N_trials x N_bands)
            % names:        frequency band names    cell (1 x 6) 
    % colors:  Contains RGB values with fields:
        % blue                                      double (1 x 3)
        % darkred                                   double (1 x 3)
    % filter_types: Cell array of strings containing fieldnames of filters
    %                                               cell (1 x 2)
    % filter_names: Cell array of strings containing titles of filters
    %                                               cell (2 x 1)
%
% OUTPUTS
    % fig:      Figure
%
% author    Miriam Kirchhoff
% project   C2B

fig = figure('units','normalized','outerposition',[0.6 0.6 0.4 0.4]); 
set(gcf,'color','w');
t = tiledlayout(2,2, 'TileSpacing','compact');
title(t, 'Relationship of MEPs and phase', 'FontName', 'Times New Roman')
lim = [-0.07 0.07];


% plot phase
for idx_filter = 1:length(filter_types)

    disp((filter_types{idx_filter}))
    nexttile
    
    %temp_size = abs(data.bandpower.(filter_types{idx_filter})(3,:))*3;
    scatter(rad2deg(data.phase.(filter_types{idx_filter})), data.mep, ...
        6, data.bandpower.(filter_types{idx_filter})(:,3)',...
        'filled', MarkerFaceAlpha= 0.6, MarkerEdgeAlpha= 0.6)
    
    % define colormap
    clear colormap_custom
    for i = 1:3
        colormap_custom(i,:) = [linspace(colors.blue(i), 1, 8), linspace(1, colors.darkred(i), 8)];
    end

    colormap(colormap_custom')
    clim([-3, 3])
    xlim(rad2deg([-pi, pi]))
    ylim([-3 2.5])
    
    grid on

    colors_temp = [[1 1 1]*0.8; colors.blue; colors.darkred];
    cutoff = 2;

    % Add fit for phase
    for i = 2:3
        if i == 1
            include_idx = logical(ones(size(data.bandpower.(filter_types{idx_filter})(:,3))));
            % disp('all')
        elseif i == 2
            include_idx = data.bandpower.(filter_types{idx_filter})(:,3) < -cutoff;
            % disp('lower')
            % disp(sum(include_idx))
        elseif i == 3
            include_idx = data.bandpower.(filter_types{idx_filter})(:,3) > cutoff;
            % disp('higher')
            % disp(sum(include_idx))
        end

        phase = double(data.phase.(filter_types{idx_filter})(include_idx));
        mep = double(data.mep(include_idx));
        T = array2table([cos(phase) sin(phase) mep]);
        T.Properties.VariableNames(:) = {'phase_cos','phase_sin', 'mep'};

        % fit linear mixed model
        fit_phase = fitlme(T, 'mep ~ phase_cos + phase_sin');

        % plot results
        x = -pi:0.01:pi;
        y = fit_phase.Coefficients.Estimate(1) + ...
            fit_phase.Coefficients.Estimate(2).* cos(x) ...
            + fit_phase.Coefficients.Estimate(3).* sin(x);
        hold on
        plot(rad2deg(x),y, 'Color',colors_temp(i,:), 'LineWidth',2)

        % write maximum location
        [max_y, max_loc] = max(y);
        max_x = rad2deg(x(max_loc));
        fprintf('max = %4.2f, max loc = %4.2f \n', max_y, max_x )
    end

    % Add oscillation
    y = cos(x)*0.5 - 2.5;
    plot(rad2deg(x),y, 'k:')

    xlabel('\mu phase [degree]')
    if idx_filter == 1
        ylabel('MEP (log, z-score)')
    else
        
        yticklabels([]);
        temp = colorbar('eastoutside');
        ylabel(temp,'\mu power','Rotation',270, 'FontSize',8)
        
    end
    title(filter_names{idx_filter})

end


% Add boxplot for peak vs trough

% plot phase
for idx_filter = 1:length(filter_types)

    disp((filter_types{idx_filter}))
    nexttile
    T = [];

    % create table for high vs low BP
    for i = 2:3
        if i == 1
            include_idx = logical(ones(size(data.bandpower.(filter_types{idx_filter})(:,3))));
            % disp('all')
        elseif i == 2
            category = "low";
            include_idx = data.bandpower.(filter_types{idx_filter})(:,3) < -cutoff;
            % disp('lower')
        elseif i == 3
            category = "high";
            include_idx = data.bandpower.(filter_types{idx_filter})(:,3) > cutoff;
            % disp('higher')
        end


        phase = double(data.phase.(filter_types{idx_filter})(include_idx));
        mep = double(data.mep(include_idx));

        peak_meps = mep(phase > deg2rad(-45) & phase < deg2rad(45));
        trough_meps = mep(phase < deg2rad(-135) | phase > deg2rad(135));

        phase_group = categorical(repelem(["peak";"trough"],[length(peak_meps),length(trough_meps)]));
        bp_group = categorical(repelem(category,length(phase_group)));
        


        temp = table([peak_meps ; trough_meps], phase_group, bp_group');
        temp.Properties.VariableNames(:) = {'mep', 'phase_group','bp_group'};
        T = [T;temp];
        
    end

    bp_group = T.("bp_group");
    p = boxchart(T.("phase_group"), T.("mep"), GroupByColor = bp_group);
    p(1).BoxFaceColor=colors.blue;
    p(2).BoxFaceColor=colors.darkred;
    %violinplot(T.("phase_group"), T.("mep"), GroupByColor = bp_group)
    

    if idx_filter == 1; ylabel('MEP (log, z-score)');
    else; yticklabels([]); end
    
    ylim([-4 3])

end

l = legend(Location="eastoutside"); title(l,'µ power');

fontsize(gcf,scale=1.4);
