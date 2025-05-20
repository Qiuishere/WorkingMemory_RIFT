function [] = plot_glm(coefficients_table, xlim_vals, xlabel_bin, variable_names_cell)
    % This function plots information from a glmm. Input required is a 
    % coefficient table (.Coefficients) from fitglme() output.

    coeffs_cell = dataset2cell(coefficients_table);
    color_vals = [1, 0, 0; 0, 0, 0; 0, 0, 1]; % red, black, blue
    sig_vals = [0.05, 0.01, 0.001]; % one, two, three stars
    
    % Prepare data
%     variable_names_cell = coeffs_cell(2:end, 1);
    betas = cell2mat(coeffs_cell(2:end, 2));
    p_vals = cell2mat(coeffs_cell(2:end, 6));
    beta_ranges = cell2mat(coeffs_cell(2:end, 7:8));

    % Update based on display order
    [~, sorted_idxs] = sort(betas, 'descend');
    variable_names_cell = variable_names_cell(sorted_idxs, 1); 
    betas = betas(sorted_idxs); p_vals = p_vals(sorted_idxs); beta_ranges = beta_ranges(sorted_idxs, :);

    % Adjust variable names
    for i = 1:length(betas)
        if variable_names_cell{i}(end-1:end) == '_1'
            variable_names(i) = string(variable_names_cell{i}(1:end-2));
        else
            variable_names(i) = string(variable_names_cell{i});
        end
    end

    % Errorbar lengths
    ebar_neg = betas - beta_ranges(:, 1); ebar_pos = beta_ranges(:, 2) - betas; 

    % Colour based on CIs
    sign_info = sum(sign(beta_ranges), 2); % -2 if neg, 0 if neutral, 2 if positive
    sign_info = ((sign_info+2)/2)+1; % -2,0,2 -> 1,2,3

    % Plot
    hold on;

    for i = 1:length(betas)
        errorbar(betas(i), i, 0, 0, ebar_neg(i), ebar_pos(i), ...
            'CapSize', 5, 'Marker', '.', 'MarkerSize', 30, 'LineWidth', 2, 'Color', color_vals(sign_info(i), :));

        % Add significance marker
        if sign_info(i) ~= 2
            num_stars = find((p_vals(i) < sig_vals) == 1, 1, 'last');
            if sign_info(i) == 1
                text(betas(i) + ebar_pos(i) + 0.001, i-0.085, repmat('*', 1, num_stars), 'FontSize', 16);
            elseif sign_info(i) == 3
                text(betas(i) - ebar_neg(i) - 0.001, i-0.075, repmat('*', 1, num_stars), 'FontSize', 16, 'HorizontalAlignment', 'right');
            end
        end
    end

    yticks([1:length(betas)]); yticklabels(variable_names); ylim([0.5 length(betas)+0.5]); xlim([xlim_vals(1), xlim_vals(2)]); ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
    xline(0, 'k--'); if xlabel_bin == 1; xlabel("Model weights", "FontName", 'Nexa-bold'); end; %ylabel("Predictors", "FontName", 'Nexa-bold'); 
    set(gca, 'TickDir', 'out');

    % Half border
    set(gca,'linewidth', 1)
    yline(length(betas)+0.5, 'LineWidth', 1); xline(xlim_vals(2), 'LineWidth', 1)
end