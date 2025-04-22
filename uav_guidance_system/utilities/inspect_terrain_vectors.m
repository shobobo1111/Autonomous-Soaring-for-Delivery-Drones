% inspect_terrain_vectors.m
% This script visualizes the terrain data and slope vectors
% to verify correct transfer from Python to MATLAB

clear; clc; close all;

% Load the terrain data exported from Python
fprintf('Loading terrain data...\n');
terrain_data = load('terrain_data.mat');

% Report basic statistics
fprintf('Data loaded.\n');
fprintf('Terrain dimensions: [%d, %d]\n', size(terrain_data.terrain));
fprintf('Elevation range: %.2f to %.2f m\n', min(min(terrain_data.terrain)), max(max(terrain_data.terrain)));
fprintf('Slope vector ranges:\n');
fprintf('  dx: %.4f to %.4f\n', min(min(terrain_data.slope_dx)), max(max(terrain_data.slope_dx)));
fprintf('  dy: %.4f to %.4f\n', min(min(terrain_data.slope_dy)), max(max(terrain_data.slope_dy)));
fprintf('  dz: %.4f to %.4f\n', min(min(terrain_data.slope_dz)), max(max(terrain_data.slope_dz)));

% Create a figure for top-down view
figure('Position', [100, 100, 1200, 800]);

% Convert lat/lon to meters for easier visualization
lat_center = mean(terrain_data.y);
lon_scale = 111000 * cos(lat_center * pi/180);  % Approx meters per degree longitude
lat_scale = 111000;  % Approx meters per degree latitude
    
% Create coordinate matrices in meters
[ny, nx] = size(terrain_data.terrain);
x_meters = (terrain_data.x - mean(terrain_data.x)) * lon_scale;
y_meters = (terrain_data.y - mean(terrain_data.y)) * lat_scale;
[X, Y] = meshgrid(x_meters, y_meters);

% Subplot 1: Terrain with contour lines
subplot(2,2,1);
contourf(X, Y, terrain_data.terrain, 20);
colormap(gca, parula);
colorbar;
title('Terrain Elevation');
xlabel('X (meters)');
ylabel('Y (meters)');
axis equal tight;

% Subplot 2: Birds-eye view with vectors
subplot(2,2,2);
% Use surf with view from above
surf(X, Y, terrain_data.terrain, 'EdgeColor', 'none');
view(2); % Top-down view
colormap(gca, parula);
hold on;

% Find local maxima (peaks) to verify vector directions
max_filter_size = ceil(min(nx, ny) / 10); % 10% of min dimension
terrain_expanded = padarray(terrain_data.terrain, [max_filter_size, max_filter_size], -inf);
[peaks_y, peaks_x] = find(imdilate(terrain_expanded, ones(2*max_filter_size+1)) == terrain_expanded);
peaks_y = peaks_y - max_filter_size;
peaks_x = peaks_x - max_filter_size;

% Only keep peaks within the valid range
valid_peaks = peaks_y >= 1 & peaks_y <= ny & peaks_x >= 1 & peaks_x <= nx;
peaks_y = peaks_y(valid_peaks);
peaks_x = peaks_x(valid_peaks);

% Keep only the top 5 peaks
if length(peaks_y) > 5
    peak_heights = zeros(length(peaks_y), 1);
    for i = 1:length(peaks_y)
        peak_heights(i) = terrain_data.terrain(peaks_y(i), peaks_x(i));
    end
    [~, idx] = sort(peak_heights, 'descend');
    peaks_y = peaks_y(idx(1:5));
    peaks_x = peaks_x(idx(1:5));
end

% Mark peaks with red circles
for i = 1:length(peaks_y)
    plot3(X(peaks_y(i), peaks_x(i)), Y(peaks_y(i), peaks_x(i)), ...
         terrain_data.terrain(peaks_y(i), peaks_x(i)) + 10, ...
         'ro', 'MarkerSize', 10, 'LineWidth', 2);
end

% Plot vectors
skip = 4; % Skip factor for visualization
quiver3(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
       terrain_data.terrain(1:skip:end, 1:skip:end), ...
       terrain_data.slope_dx(1:skip:end, 1:skip:end), ...
       terrain_data.slope_dy(1:skip:end, 1:skip:end), ...
       zeros(size(terrain_data.terrain(1:skip:end, 1:skip:end))), ... 
       0.5, 'b', 'LineWidth', 1);
title('Terrain with Horizontal Slope Vectors');
xlabel('X (meters)');
ylabel('Y (meters)');

% Subplot 3: 3D perspective with vectors
subplot(2,2,3:4);
surf(X, Y, terrain_data.terrain, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold on;

% Plot full 3D vectors
quiver3(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
       terrain_data.terrain(1:skip:end, 1:skip:end), ...
       terrain_data.slope_dx(1:skip:end, 1:skip:end), ...
       terrain_data.slope_dy(1:skip:end, 1:skip:end), ...
       terrain_data.slope_dz(1:skip:end, 1:skip:end), ...
       0.5, 'r', 'LineWidth', 1);

% Mark peaks
for i = 1:length(peaks_y)
    plot3(X(peaks_y(i), peaks_x(i)), Y(peaks_y(i), peaks_x(i)), ...
         terrain_data.terrain(peaks_y(i), peaks_x(i)) + 10, ...
         'ro', 'MarkerSize', 10, 'LineWidth', 2);
end

title('3D View with Slope Vectors');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Elevation (m)');
view(30, 30);
grid on;

% Display instructions
fprintf('\nInterpreting the plot:\n');
fprintf('- Red circles mark the highest peaks\n');
fprintf('- Blue arrows in top-down view show horizontal components (dx, dy)\n');
fprintf('- Red arrows in 3D view show full slope vectors (dx, dy, dz)\n');
fprintf('- For correct upslope vectors, arrows should point toward peaks\n');