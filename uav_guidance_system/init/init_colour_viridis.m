function [custom_cmap, named_colors] = init_colour_viridis()
    % INIT_COLOUR_VIRIDIS Initialses Viridis colormap for plotting
    %
    % Returns:
    %   custom_cmap: 100-element Viridis colormap for path coloring
    %   named_colors: Structure with named colors from the palette
    
    % Define viridis colormap colors with meaningful names
    % These are key colors from the Viridis colormap
    color_names = {'dark_purple', 'purple', 'blue', 'teal', 'green', 'lime', 'yellow'};
    
    % Viridis hex colors at different points in the spectrum
    custom_colors = {
        '#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#bddf26', '#fde725'
    };
    
    % Convert hex to RGB
    custom_rgb = zeros(length(custom_colors), 3);
    for i = 1:length(custom_colors)
        hex = custom_colors{i}(2:end);
        custom_rgb(i, 1) = hex2dec(hex(1:2))/255; % Red
        custom_rgb(i, 2) = hex2dec(hex(3:4))/255; % Green
        custom_rgb(i, 3) = hex2dec(hex(5:6))/255; % Blue
        
        % Store in named_colors structure
        named_colors.(color_names{i}) = custom_rgb(i,:);
    end
    
    % Interpolate to create a 100-step colormap
    custom_cmap = interp1(linspace(0, 1, size(custom_rgb, 1)), custom_rgb, linspace(0, 1, 100));
    
    % Add additional named colors for specific uses
    named_colors.light = [0.95, 0.95, 0.95];   % Almost white
    named_colors.dark = [0.2, 0.2, 0.2];       % Dark gray for text
    named_colors.highlight = [0.85, 0.37, 0];  % Orange highlight color
end