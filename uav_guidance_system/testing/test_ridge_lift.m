% test_ridge_lift.m
clear; clc;

% Load terrain data
terrain_data = load('terrain_data.mat');

% Print diagnostic information
fprintf('\n===== TERRAIN AND SLOPE DIAGNOSTICS =====\n');
fprintf('Terrain dimensions: [%d, %d]\n', size(terrain_data.terrain));
fprintf('Terrain elevation range: %.2f to %.2f m\n', min(min(terrain_data.terrain)), max(max(terrain_data.terrain)));

% Convert lat/lon to meters for consistent visualization
lat_center = mean(terrain_data.y);
lon_scale = 111000 * cos(lat_center * pi/180);  % Meters per degree longitude
lat_scale = 111000;  % Meters per degree latitude

% Create coordinate vectors in meters
x_meters = (terrain_data.x - mean(terrain_data.x)) * lon_scale;
y_meters = (terrain_data.y - mean(terrain_data.y)) * lat_scale;

% Create wind field with meshgrid for consistency
fprintf('\n===== INITIALISING WIND FIELD =====\n');
x_range = [min(x_meters), max(x_meters)];
y_range = [min(y_meters), max(y_meters)];
z_range = [0, 500];  % 0-500m above ground level
grid_spacing = (max(x_meters) - min(x_meters)) / (length(x_meters) - 1);
fprintf('Grid spacing: %.2f m\n', grid_spacing);

% Initialize the wind field
wind_field = init_wind_field(x_range, y_range, z_range, grid_spacing);

% Set up basic wind parameters
wind_speed = 2; % m/s
wind_direction = 100; % degrees (East)
fprintf('Wind: %.1f m/s from %.0f degrees\n', wind_speed, wind_direction);

% Create terrain visualization before adding ridge lift
fprintf('\n===== VISUALIZING TERRAIN AND SLOPE VECTORS =====\n');
figure('Position', [100, 500, 900, 400]);

% First: plot terrain with pure slope vectors for verification
subplot(1,2,1);
% Create meshgrid coordinates for terrain
[X_terrain, Y_terrain] = meshgrid(x_meters, y_meters);
surf(X_terrain, Y_terrain, terrain_data.terrain, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
colormap(gca, 'parula');
hold on;

% Plot slope vectors on terrain
skip = 4; % Spacing for vectors
vector_scale = 20; % Scale for better visibility

% Use the terrain coordinates directly for vector starting points
quiver3(X_terrain(1:skip:end, 1:skip:end), ...
       Y_terrain(1:skip:end, 1:skip:end), ...
       terrain_data.terrain(1:skip:end, 1:skip:end), ...
       terrain_data.slope_dx(1:skip:end, 1:skip:end) * vector_scale, ...
       terrain_data.slope_dy(1:skip:end, 1:skip:end) * vector_scale, ...
       terrain_data.slope_dz(1:skip:end, 1:skip:end) * vector_scale, ...
       0, 'r', 'LineWidth', 1.5);

title('Terrain with Slope Vectors');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevation (m)');
view(30, 30);

% Apply ridge lift to wind field
fprintf('\n===== APPLYING RIDGE LIFT =====\n');
wind_field = create_ridge_lift(terrain_data, wind_speed, wind_direction, wind_field);

% Visualize wind field with ridge lift
subplot(1,2,2);
surf(X_terrain, Y_terrain, terrain_data.terrain, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;

% Add wind vectors at two different heights
level1 = 1;  % Ground level
level2 = 5;  % Higher altitude (depends on your z spacing)

% Sample the wind field at these levels
% Make sure to use the smaller of terrain and wind field dimensions
[ny, nx] = size(terrain_data.terrain);
[wy, wx, wz] = size(wind_field.X);
use_ny = min(ny, wy);
use_nx = min(nx, wx);

% Scale for wind vectors
wind_scale = 200;


% Plot wind vectors at ground level
quiver3(wind_field.X(1:skip:use_ny, 1:skip:use_nx, level1), ...
       wind_field.Y(1:skip:use_ny, 1:skip:use_nx, level1), ...
       wind_field.Z(1:skip:use_ny, 1:skip:use_nx, level1), ...
       wind_field.Wx(1:skip:use_ny, 1:skip:use_nx, level1) * wind_scale, ...
       wind_field.Wy(1:skip:use_ny, 1:skip:use_nx, level1) * wind_scale, ...
       wind_field.U(1:skip:use_ny, 1:skip:use_nx, level1) * wind_scale, ...
       0, 'b', 'LineWidth', 1.5);

% Plot wind vectors at higher level 
quiver3(wind_field.X(1:skip:use_ny, 1:skip:use_nx, level2), ...
       wind_field.Y(1:skip:use_ny, 1:skip:use_nx, level2), ...
       wind_field.Z(1:skip:use_ny, 1:skip:use_nx, level2), ...
       wind_field.Wx(1:skip:use_ny, 1:skip:use_nx, level2) * wind_scale, ...
       wind_field.Wy(1:skip:use_ny, 1:skip:use_nx, level2) * wind_scale, ...
       wind_field.U(1:skip:use_ny, 1:skip:use_nx, level2) * wind_scale, ...
       0, 'g', 'LineWidth', 1.5);

title('Ridge Lift Wind Field');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevation (m)');
view(30, 30);

% Add legend
h1 = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 1.5);
h2 = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 1.5);
h3 = plot3(NaN, NaN, NaN, 'g-', 'LineWidth', 1.5);
legend([h1, h2, h3], {'Terrain Slope', 'Wind (Ground)', 'Wind (Higher)'}, 'Location', 'northeast');

% Add cross-section view for better understanding
figure('Position', [100, 100, 900, 300]);
% Take a slice through the middle of y-dimension
mid_y = round(use_ny/2);
x_indices = 1:use_nx;

% Plot terrain profile
plot(x_meters(x_indices), terrain_data.terrain(mid_y, x_indices), 'k-', 'LineWidth', 2);
hold on;

% Add slope vectors to terrain
skip_profile = 8;  % Increase skip for clarity
x_pos = x_meters(1:skip_profile:use_nx);
y_pos = terrain_data.terrain(mid_y, 1:skip_profile:use_nx);
sx = terrain_data.slope_dx(mid_y, 1:skip_profile:use_nx);
sz = terrain_data.slope_dz(mid_y, 1:skip_profile:use_nx);
quiver(x_pos, y_pos, sx*vector_scale, sz*vector_scale, 0, 'r', 'LineWidth', 1.5);

% Add wind vectors at two heights
height1 = 50;  % Meters above terrain
height2 = 150; % Meters above terrain

% Find the closest z-levels to these heights
[~, z_idx1] = min(abs(wind_field.z - height1));
[~, z_idx2] = min(abs(wind_field.z - height2));

% Extract wind data at these heights
y1 = y_pos + height1;
u1 = wind_field.U(mid_y, 1:skip_profile:use_nx, z_idx1);
wx1 = wind_field.Wx(mid_y, 1:skip_profile:use_nx, z_idx1);

y2 = y_pos + height2;
u2 = wind_field.U(mid_y, 1:skip_profile:use_nx, z_idx2);
wx2 = wind_field.Wx(mid_y, 1:skip_profile:use_nx, z_idx2);

% Plot wind vectors
vertical_scale = 5;  % Adjust for vertical wind component visibility
horizontal_scale = 5;  % Adjust for horizontal wind component visibility
quiver(x_pos, y1, wx1*horizontal_scale, u1*vertical_scale, 0, 'b', 'LineWidth', 1.5);
quiver(x_pos, y2, wx2*horizontal_scale, u2*vertical_scale, 0, 'g', 'LineWidth', 1.5);

% Enhance the visualization
title('Terrain Cross-Section with Slope and Wind Vectors');
xlabel('X distance (m)');
ylabel('Elevation (m)');
grid on;
legend({'Terrain', 'Slope Vectors', sprintf('Wind at %dm', height1), sprintf('Wind at %dm', height2)}, 'Location', 'northeast');

% Add some headroom to the plot
ylim([min(min(terrain_data.terrain))-50, max(max(terrain_data.terrain))+height2+100]);

fprintf('\nVisualization complete. Check the figures to verify slope vectors and wind field.\n');



% % test_ridge_lift.m
% clear; clc;
% 
% % Load terrain data
% terrain_data = load('terrain_data.mat');
% 
% 
% % % Print diagnostic information about slope vectors
% % fprintf('\n===== TERRAIN AND SLOPE DIAGNOSTICS =====\n');
% % fprintf('Terrain dimensions: [%d, %d]\n', size(terrain_data.terrain));
% % fprintf('Terrain elevation range: %.2f to %.2f m\n', min(min(terrain_data.terrain)), max(max(terrain_data.terrain)));
% % 
% % % Check if we have the new slope format
% % if isfield(terrain_data, 'slope_dx')
% %     fprintf('Using new slope vector format\n');
% %     % Basic statistics on slope vectors
% %     fprintf('Slope X component range: %.4f to %.4f\n', min(min(terrain_data.slope_dx)), max(max(terrain_data.slope_dx)));
% %     fprintf('Slope Y component range: %.4f to %.4f\n', min(min(terrain_data.slope_dy)), max(max(terrain_data.slope_dy)));
% %     fprintf('Slope Z component range: %.4f to %.4f\n', min(min(terrain_data.slope_dz)), max(max(terrain_data.slope_dz)));
% % 
% %     % Check for NaN or Inf values
% %     fprintf('NaN values in slope_dx: %d\n', sum(sum(isnan(terrain_data.slope_dx))));
% %     fprintf('Inf values in slope_dx: %d\n', sum(sum(isinf(terrain_data.slope_dx))));
% % 
% %     % Check vector normalization
% %     mag_check = sqrt(terrain_data.slope_dx.^2 + terrain_data.slope_dy.^2 + terrain_data.slope_dz.^2);
% %     fprintf('Average vector magnitude: %.4f (should be close to 1)\n', mean(mag_check(:)));
% % else
% %     fprintf('WARNING: Using old slope vector format - results may be incorrect\n');
% % end
% 
% % Create wind field
% wind_speed = 2; % m/s
% wind_direction = 90; % degrees (North)
% wind_field = create_ridge_lift(terrain_data, wind_speed, wind_direction);
% 
% 
% % % Check correlation between slope and vertical wind at surface
% % level1 = 1;
% % if isfield(terrain_data, 'slope_dx')
% %     % Get slope and vertical wind at surface
% %     surface_wind = wind_field.U(:,:,level1);
% %     slope_x = terrain_data.slope_dx;
% %     slope_y = terrain_data.slope_dy;
% % 
% %     % Calculate alignment with wind direction
% %     wind_dir_rad = (270 - wind_direction) * pi/180;
% %     wind_unit = [cos(wind_dir_rad), sin(wind_dir_rad)];
% % 
% %     % Compute dot product for each point
% %     alignment = zeros(size(slope_x));
% %     for i = 1:size(slope_x, 1)
% %         for j = 1:size(slope_x, 2)
% %             slope_vec = [slope_x(i,j), slope_y(i,j)];
% %             slope_mag = norm(slope_vec);
% %             if slope_mag > 0.001
% %                 slope_unit = slope_vec / slope_mag;
% %                 alignment(i,j) = dot(wind_unit, slope_unit);
% %             end
% %         end
% %     end
% % 
% %     fprintf('Alignment with wind: average %.4f, range %.4f to %.4f\n', mean(alignment(:)), min(min(alignment)), max(max(alignment)));
% %     fprintf('Average wind deflection at surface: %.4f m/s\n', mean(surface_wind(:)));
% % 
% %     % Check correlation
% %     [r, p] = corrcoef(alignment(:), surface_wind(:));
% %     fprintf('Correlation between alignment and vertical wind: r=%.4f, p=%.4f\n', r(1,2), p(1,2));
% % end
% 
% % Create figure
% figure('Position', [100, 100, 1200, 600]);
% 
% % Plot terrain surface
% subplot(1,2,1);
% surf(wind_field.X(:,:,1), wind_field.Y(:,:,1), terrain_data.terrain, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.8);
% colormap(gca, 'parula');
% title('Terrain Surface');
% xlabel('X distance (m)');
% ylabel('Y distance (m)');
% zlabel('Elevation (m)');
% view(45, 30);
% 
% % Wind field plot with slope vectors
% subplot(1,2,2);
% 
% % Plot terrain surface with increased transparency
% surf(wind_field.X(:,:,1), wind_field.Y(:,:,1), terrain_data.terrain, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.2);
% hold on;
% 
% % Spacing for vectors (horizontal)
% n = 2; 
% 
% % Select just two vertical levels to plot
% level1 = 1;  % First level (right at terrain surface)
% level2 = 10; % Second level (higher above terrain)
% 
% % Plot the first level of vectors (at the slope)
% x = wind_field.X(1:n:end, 1:n:end, level1);
% y = wind_field.Y(1:n:end, 1:n:end, level1);
% z = wind_field.Z(1:n:end, 1:n:end, level1);
% wx = wind_field.Wx(1:n:end, 1:n:end, level1);
% wy = wind_field.Wy(1:n:end, 1:n:end, level1);
% u = wind_field.U(1:n:end, 1:n:end, level1);
% 
% % Scale vectors for visibility
% scale = 20;
% quiver3(x, y, z, wx*scale, wy*scale, u*scale, 0, ...
%     'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
% 
% % Plot the second level of vectors (higher above terrain)
% x = wind_field.X(1:n:end, 1:n:end, level2);
% y = wind_field.Y(1:n:end, 1:n:end, level2);
% z = wind_field.Z(1:n:end, 1:n:end, level2);
% wx = wind_field.Wx(1:n:end, 1:n:end, level2);
% wy = wind_field.Wy(1:n:end, 1:n:end, level2);
% u = wind_field.U(1:n:end, 1:n:end, level2);
% 
% % Scale vectors for visibility
% quiver3(x, y, z, wx*scale, wy*scale, u*scale, 0, ...
%     'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
% 
% % Plot slope vectors if using new format
% if isfield(terrain_data, 'slope_dx')
%     x_slope = wind_field.X(1:n:end, 1:n:end, level1);
%     y_slope = wind_field.Y(1:n:end, 1:n:end, level1);
%     z_slope = wind_field.Z(1:n:end, 1:n:end, level1);
% 
%     % Get slope data
%     sx = terrain_data.slope_dx(1:n:end, 1:n:end);
%     sy = terrain_data.slope_dy(1:n:end, 1:n:end);
%     sz = terrain_data.slope_dz(1:n:end, 1:n:end);
% 
%     % Plot terrain slope vectors
%     slope_scale = 20; 
%     quiver3(x_slope, y_slope, z_slope, sx*slope_scale, sy*slope_scale, sz*slope_scale, 0, ...
%         'Color', [0, 0.7, 0], 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
% 
%     fprintf('\nSlope vectors plotted in green\n');
% end
% 
% % Adjust aspect ratio and view
% daspect([1 1 0.2]);
% view(30, 30);
% title('Wind Field with Terrain Slope Vectors');
% xlabel('X distance (m)');
% ylabel('Y distance (m)');
% zlabel('Height (m)');
% grid on;
% 
% % Set reasonable axis limits
% xlim([min(wind_field.x) max(wind_field.x)]);
% ylim([min(wind_field.y) max(wind_field.y)]);
% zlim([min(min(terrain_data.terrain)) max(max(terrain_data.terrain))+500]);
% 
% % Add legend
% h1 = plot3(NaN, NaN, NaN, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
% h2 = plot3(NaN, NaN, NaN, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
% h3 = plot3(NaN, NaN, NaN, 'Color', [0, 0.7, 0], 'LineWidth', 2);
% legend([h1, h2, h3], {'Wind at Surface', 'Wind at Height', 'Terrain Slope'}, 'Location', 'northeast');
% 
% % Create a cross-section view for better analysis
% figure('Position', [150, 150, 900, 400]);
% mid_y = round(size(terrain_data.terrain, 1)/2);
% x_indices = 1:size(terrain_data.terrain, 2);
% 
% % Plot the terrain profile
% plot(wind_field.x, terrain_data.terrain(mid_y,:), 'k-', 'LineWidth', 2);
% hold on;
% 
% % Sample slope vectors along profile
% slope_sample = 4;
% x_pos = wind_field.x(1:slope_sample:end);
% y_pos = terrain_data.terrain(mid_y,1:slope_sample:end);
% 
% if isfield(terrain_data, 'slope_dx')
%     sx = terrain_data.slope_dx(mid_y,1:slope_sample:end);
%     sy = terrain_data.slope_dz(mid_y,1:slope_sample:end);  % Use dz as y for 2D plot
% 
%     % Scale for display
%     slope_scale_2d = 100;
%     quiver(x_pos, y_pos, sx*slope_scale_2d, sy*slope_scale_2d, 0, 'g', 'LineWidth', 1.5);
% 
%     % Add wind vectors at two heights
%     height1 = 10;  % Meters above terrain
%     height2 = 10; % Meters above terrain
% 
%     % Get wind at first height
%     y1 = y_pos + height1;
%     u1 = squeeze(wind_field.U(mid_y,1:slope_sample:end,2));  % Level 2 approximates 20m
% 
%     % Get wind at second height
%     y2 = y_pos + height2;
%     u2 = squeeze(wind_field.U(mid_y,1:slope_sample:end,6));  % Level 6 approximates 100m
% 
%     % Plot wind vectors (vertical component only for clarity)
%     wind_scale = 20;
%     quiver(x_pos, y1, zeros(size(x_pos)), u1*wind_scale, 0, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
%     quiver(x_pos, y2, zeros(size(x_pos)), u2*wind_scale, 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
% 
%     % Add legend
%     h1 = plot(NaN, NaN, 'g-', 'LineWidth', 2);
%     h2 = plot(NaN, NaN, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.002);
%     h3 = plot(NaN, NaN, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.002);
%     legend([h1, h2, h3], {'Terrain Slope', sprintf('Wind at %dm', height1), sprintf('Wind at %dm', height2)}, 'Location', 'northeast');
% end
% 
% title('Terrain Cross-Section with Slope and Wind Vectors');
% xlabel('X distance (m)');
% ylabel('Elevation (m)');
% grid on;
% 
% fprintf('\nVisualization complete. Check the figures for terrain, slope, and wind patterns.\n');