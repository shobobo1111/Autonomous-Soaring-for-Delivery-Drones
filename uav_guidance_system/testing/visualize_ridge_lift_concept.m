function visualize_ridge_lift_concept()
    % VISUALIZE_RIDGE_LIFT_CONCEPT Demonstrates how wind interacts with terrain slope
    % Creates a 3D visualization showing wind vector transformation over terrain
    
    % Create a 3D terrain surface (Gaussian hill)
    [X, Y] = meshgrid(linspace(-100, 100, 30));
    Z = 30 * exp(-((X/40).^2 + (Y/40).^2));
    
    % Create figure with 3D view
    figure('Position', [100, 100, 900, 700]);
    
    % Plot the terrain surface
    surf(X, Y, Z, 'FaceAlpha', 0.8, 'EdgeAlpha', 0.2, 'FaceColor', [0.8, 0.9, 0.8]);
    hold on;
    
    % Calculate terrain gradients
    [dZdx, dZdy] = gradient(Z, X(1,2)-X(1,1), Y(2,1)-Y(1,1));
    
    % Set wind parameters
    wind_speed = 10;  % m/s
    wind_dir = 270;   % Wind from the west
    wind_dir_rad = (270 - wind_dir) * pi/180;
    wind_vector = [cos(wind_dir_rad), sin(wind_dir_rad), 0] * wind_speed;  % [x,y,z]
    wind_unit_horiz = [cos(wind_dir_rad), sin(wind_dir_rad)];  % Unit vector (horizontal only)
    
    % Height-based decay parameters
    decay_height = 40;  % Height where effect reduces to ~37%
    vertical_strength = 1.0;  % Multiplier for vertical component
    
    % Sample points for vector visualization
    step = 4;  % Sample every Nth point
    sample_heights = [10, 20, 30, 50];  % Heights above terrain
    
    % Plot vectors at sample points
    for i = 1:step:size(X, 1)
        for j = 1:step:size(X, 2)
            x_pos = X(i,j);
            y_pos = Y(i,j);
            terrain_height = Z(i,j);
            
            % Get terrain slope at this point
            slope_x = dZdx(i,j);
            slope_y = dZdy(i,j);
            
            % Create upslope vector (horizontal with vertical component based on slope)
            % This represents the direction a ball would roll uphill
            upslope_vector = [slope_x, slope_y];  % Horizontal components only
            upslope_mag = norm(upslope_vector);
            
            if upslope_mag > 0.001  % Only if there's a meaningful slope
                % Normalize upslope vector
                upslope_unit = upslope_vector / upslope_mag;
                
                % Create 3D display vector with small vertical component
                upslope_display = [upslope_unit(1), upslope_unit(2), 0.2] * 10;
                
                % Plot upslope vector at surface
                quiver3(x_pos, y_pos, terrain_height, ...
                       upslope_display(1), upslope_display(2), upslope_display(3), ...
                       'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
                
                % Calculate and visualize wind vectors at different heights
                for h = 1:length(sample_heights)
                    height = sample_heights(h);
                    z_pos = terrain_height + height;
                    
                    % Calculate decay factor based on height
                    decay_factor = exp(-height/decay_height);
                    
                    % Calculate dot product between wind and upslope vector
                    % Both wind_unit_horiz and upslope_unit are 2D vectors
                    dot_product = dot(wind_unit_horiz, upslope_unit);
                    
                    % Calculate vertical deflection
                    vertical_component = dot_product * vertical_strength * wind_speed * decay_factor;
                    
                    % Create resultant vector
                    resultant = [wind_vector(1), wind_vector(2), vertical_component];
                    
                    % Scale for visualization
                    scale = 7;
                    
                    % Determine color based on height (fade from red to light blue)
                    color = [1-h/length(sample_heights), 0.4, 0.4+h/length(sample_heights)*0.6];
                    
                    % Plot resultant wind vector
                    quiver3(x_pos, y_pos, z_pos, ...
                           resultant(1)/wind_speed*scale, resultant(2)/wind_speed*scale, resultant(3)/wind_speed*scale, ...
                           'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
                end
            end
        end
    end
    
    % Add a legend arrow for ambient wind
    quiver3(80, 0, 70, -20, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.8);
    text(80, 0, 75, 'Ambient Wind', 'Color', 'blue', 'FontWeight', 'bold');
    
    % Add upslope vector label
    text(30, 30, 45, 'Upslope Direction', 'Color', 'green', 'FontWeight', 'bold');
    
    % Add resultant vector label
    text(0, -50, 60, 'Resultant Wind', 'Color', [0.8, 0.3, 0.5], 'FontWeight', 'bold');
    
    % Set axis properties
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Height (meters)');
    title('Ridge Lift: 3D Vector Transformation', 'FontSize', 16);
    grid on;
    
    % Set optimal viewing angle
    view(140, 30);
    
    % Create the detailed explanation figure
    figure('Position', [200, 200, 900, 600]);
    
    % Create a more detailed terrain profile for 2D demonstration
    x_detailed = linspace(-100, 100, 201);
    z_detailed = 30 * exp(-(x_detailed/30).^2);  % Steeper Gaussian hill
    dx_detailed = gradient(z_detailed, x_detailed);
    
    % Plot the terrain
    plot(x_detailed, z_detailed, 'k', 'LineWidth', 2);
    hold on;
    
    % Select one point for demonstration
    x_demo = -20;  % Point on the windward slope
    [~, idx_demo] = min(abs(x_detailed - x_demo));
    terrain_height_demo = z_detailed(idx_demo);
    terrain_slope_demo = dx_detailed(idx_demo);
    
    % Mark the demonstration point
    plot(x_demo, terrain_height_demo, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Calculate upslope direction at demo point
    upslope_vector = [1, terrain_slope_demo];  
    upslope_mag = norm(upslope_vector);
    upslope_unit = upslope_vector / upslope_mag;
    
    % Show the upslope vector
    quiver(x_demo, terrain_height_demo, upslope_unit(1)*10, upslope_unit(2)*10, 0, 'g', 'LineWidth', 2);
    text(x_demo+12, terrain_height_demo+2, 'Upslope Direction', 'Color', 'g', 'FontWeight', 'bold');
    
    % Show the wind vector
    arrow_start = [x_demo-20, terrain_height_demo+20];
    quiver(arrow_start(1), arrow_start(2), wind_vector(1)*1.5, 0, 0, 'b', 'LineWidth', 2);
    text(arrow_start(1)-5, arrow_start(2)+5, 'Wind Direction', 'Color', 'b', 'FontWeight', 'bold');
    
    % Calculate and show the resulting vectors at different heights
    heights = [5, 15, 30, 60];
    colors = copper(length(heights));
    
    for h = 1:length(heights)
        height = heights(h);
        z_pos = terrain_height_demo + height;
        
        % Calculate decay factor
        decay_factor = exp(-height/decay_height);
        
        % Calculate dot product between wind and upslope (2D only)
        dot_product = dot([1, 0], upslope_unit);  
        
        % Calculate vertical multiplier
        vertical_multiplier = dot_product * vertical_strength;
        
        % Create resultant vector
        resultant = [wind_vector(1), vertical_multiplier * decay_factor * wind_speed];
        
        % Draw the resultant vector
        quiver(x_demo, z_pos, resultant(1)/wind_speed*15, resultant(2)/wind_speed*15, 0, ...
               'Color', colors(h,:), 'LineWidth', 2);
        
        % Add text label with height and vertical component
        text(x_demo+15, z_pos, sprintf('Height: %dm, V_z: %.2f m/s', height, resultant(2)), ...
             'Color', colors(h,:), 'FontWeight', 'bold');
    end
    
    % Add explanation text
    textbox = annotation('textbox', [0.02, 0.02, 0.4, 0.25], 'String', {...
        'Ridge Lift Calculation:', ...
        '', ...
        '1. Calculate angle between wind and upslope direction', ...
        '2. Higher dot product = stronger upward deflection', ...
        '3. Effect decreases exponentially with height', ...
        sprintf('4. At height=%.0fm: deflection reduced to %.0f%%', decay_height, 100/exp(1))}, ...
        'FitBoxToText', 'on', 'BackgroundColor', [1, 1, 0.9], ...
        'EdgeColor', [0.7, 0.7, 0.7], 'FontSize', 11);
    
    % Set axis properties
    axis([-100, 100, 0, 80]);
    xlabel('Distance (m)');
    ylabel('Height (m)');
    title('Ridge Lift Cross-Section - Westerly Wind', 'FontSize', 14);
    grid on;
end