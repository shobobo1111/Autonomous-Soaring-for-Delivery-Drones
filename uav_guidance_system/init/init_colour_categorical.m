function [custom_cmap, named_colors] = init_colour_categorical()
    % INIT_COLOUR_CATEGORICAL Initialize categorical color scheme for mode visualization
    %
    % Returns:
    %   custom_cmap: A 100-value interpolated colormap for continuous plotting
    %   named_colors: Structure with named colors for specific plot elements
    
    % Define categorical colors for operating modes (as char arrays, not strings)
    % categorical_colors = {
    %     'rgb(69,83,194)',    % Mode 1: Loiter - Blue
    %     'rgb(70,189,179)',   % Mode 2: Urgent - Teal
    %     'rgb(0,121,97)',     % Mode 3: Long Distance - Green
    %     'rgb(12,168,46)',    % Mode 4: Balanced - Light Green
    %     'rgb(49,99,135)',    % Mode 5: Energy Saving - Navy Blue
    %     'rgb(131,172,243)'   % Mode 6: Regeneration - Light Blue
    % };


    categorical_colors = {
    'rgb(82,91,118)',    % Mode 1: Loiter --purple
    'rgb(213,62,79)',   % Mode 2: Urgent - red
    'rgb(254,224,139)',     % Mode 3: Long Distance --- no lonnger used
    'rgb(252,141,89)',    % Mode 4: (balanced) thermal climb - orange
    'rgb(50,136,189)',    % Mode 5: Energy Saving -- blue
    'rgb(153,213,148)'   % Mode 6: Regeneration  --- green
    };


    % Convert RGB strings to numeric RGB values
    mode_colors = zeros(length(categorical_colors), 3);
    for i = 1:length(categorical_colors)
        % Extract RGB values from the string using regular expressions
        rgb_str = categorical_colors{i};
        rgb_parts = regexp(rgb_str, 'rgb\((\d+),(\d+),(\d+)\)', 'tokens');
        
        if ~isempty(rgb_parts)
            r = str2double(rgb_parts{1}{1});
            g = str2double(rgb_parts{1}{2});
            b = str2double(rgb_parts{1}{3});
            mode_colors(i,:) = [r, g, b] / 255; % Normalize to 0-1 range
        else
            % Default fallback if parsing fails
            mode_colors(i,:) = [0.5, 0.5, 0.5];
            warning('Failed to parse color: %s', rgb_str);
        end
    end
    
    % Create a structure with named colors for easy reference
    named_colors = struct();
    
    % Add mode colors to the structure
    mode_names = {'loiter', 'urgent', 'long_distance', 'balanced', 'energy_saving', 'regeneration'};
    for i = 1:length(mode_names)
        named_colors.(mode_names{i}) = mode_colors(i,:);
    end
    
    % Add some additional colors for other plot elements
    named_colors.yellow = [0.9290, 0.6940, 0.1250];  % For velocity difference
    named_colors.orange = [0.8500, 0.3250, 0.0980];  % For battery level
    named_colors.green = [0.4660, 0.6740, 0.1880];   % For altitude
    named_colors.blue = [0.3010, 0.7450, 0.9330];    % General blue
    named_colors.red = [0.6350, 0.0780, 0.1840];     % For warnings/highlights
    
    % Create a 100-value custom colormap for continuous plotting
    custom_cmap = interp1(linspace(0, 1, size(mode_colors, 1)), mode_colors, linspace(0, 1, 100));
end