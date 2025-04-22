function [custom_cmap, named_colors] = init_colour()
    % Initialise custom colors for plotting
    %
    % Returns:
    %   custom_cmap: 100-element colormap for path coloring
    %   blue_color: RGB triplet for ground speed line
    %   red_color: RGB triplet for airspeed line
    %   named_colors: Structure with named colors from the palette

    % Define spectral colormap with color names
    color_names = {'dark_red', 'red', 'orange', 'light_orange', 'yellow', ...
                  'light_green', 'mint_green', 'teal', 'blue', 'purple'};
    
    custom_colors = {
        '#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', ...
        '#e6f598', '#abdda4', '#66c2a5', '#3288bd', '#5e4fa2'
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
end