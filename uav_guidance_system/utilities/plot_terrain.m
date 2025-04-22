% % plot_terrain_standalone.m
% % Self-contained script to visualize terrain data with a green-to-brown color scheme
% % Just run this script - no function call required
% 
% % Specify the path to your terrain data file
% terrain_data_filename = 'terrain_data.mat';
% 
% % Load terrain data
% try
%     terrain_data = load(terrain_data_filename);
%     fprintf('Terrain data loaded successfully from: %s\n', terrain_data_filename);
% catch err
%     error('Error loading terrain data: %s', err.message);
% end
% 
% % Set up global plotting settings for consistency
% set(0, 'DefaultAxesFontSize', 12);
% set(0, 'DefaultFigureColor', 'w');
% set(0, 'DefaultTextInterpreter', 'tex');
% set(0, 'DefaultAxesFontName', 'Times New Roman');
% set(0, 'DefaultTextFontName', 'Times New Roman');
% 
% % Create figure with thesis-appropriate dimensions (half width)
% fig_handle = figure;
% set(fig_handle, 'Units', 'centimeters');
% set(fig_handle, 'Position', [0 0 8.5 5]); % Half of the 17cm width
% 
% % Get coordinate conversion factors
% lat_center = mean(terrain_data.y);
% lon_scale = 111000 * cos(lat_center * pi/180);  % meters per degree longitude
% lat_scale = 111000;  % meters per degree latitude
% 
% % Convert to meters
% x_meters = (terrain_data.x - mean(terrain_data.x)) * lon_scale;
% y_meters = (terrain_data.y - mean(terrain_data.y)) * lat_scale;
% 
% % Create terrain mesh
% [X_terrain, Y_terrain] = meshgrid(x_meters, y_meters);
% 
% % Plot terrain surface with green to brown colormap
% terrain_handle = surf(X_terrain, Y_terrain, terrain_data.terrain, 'EdgeColor', 'none', 'FaceAlpha', 1);
% 
% % Create custom green to brown colormap
% green_to_brown = [
%     0.2, 0.6, 0.2;   % Dark green (low)
%     0.3, 0.7, 0.3;   % Medium green
%     0.5, 0.8, 0.3;   % Light green
%     0.7, 0.9, 0.4;   % Yellow-green
%     0.8, 0.7, 0.3;   % Light brown
%     0.6, 0.5, 0.2;   % Medium brown
%     0.5, 0.3, 0.1    % Dark brown (high)
% ];
% 
% % Apply custom colormap
% colormap(gca, green_to_brown);
% 
% % Add colorbar
% c = colorbar;
% ylabel(c, 'Elevation (m)', 'FontWeight', 'bold');
% 
% % Enhance lighting
% lighting gouraud
% material dull
% light('Position', [0.1 0.4 0.4], 'Style', 'infinite');
% 
% % Format axes and labels
% xlabel('East (m)', 'FontWeight', 'bold');
% ylabel('North (m)', 'FontWeight', 'bold');
% zlabel('Elevation (m)', 'FontWeight', 'bold');
% 
% % Set position of axes in the figure for consistency
% set(gca, 'Units', 'centimeters');
% set(gca, 'Position', [1.5 1.5 6 10]); % Adjusted for half-width
% 
% % Set view and grid
% grid off;
% view(20, 20);
% 
% % Set axis limits and scaling
% axis tight;
% axis equal;
% 
% % Box on for cleaner appearance
% box on;
% 
% % Add title
% title('Wicklow Mountains Terrain', 'FontWeight', 'bold');
% 
% % Print completion message
% fprintf('Terrain visualization complete\n');



% plot_terrain_improved.m
% Self-contained script to visualize terrain data with improved continuous color scheme
% Just run this script - no function call required

%% USER ADJUSTABLE PARAMETERS %%
% Specify the path to your terrain data file
terrain_data_filename = 'terrain_data.mat';

% Figure dimensions in centimeters
fig_width = 8.5;  % Width in cm (half of 17cm)
fig_height = 5;  % Height in cm

% Plot padding within figure (in cm)
left_padding = 1.5;
bottom_padding = 1.5;
right_padding = 1.0;
top_padding = 1.0;

% View angle
view_azimuth = 20;  % Horizontal rotation
view_elevation = 20;  % Vertical elevation angle

%% SCRIPT EXECUTION %%
% Load terrain data
try
    terrain_data = load(terrain_data_filename);
    fprintf('Terrain data loaded successfully from: %s\n', terrain_data_filename);
catch err
    error('Error loading terrain data: %s', err.message);
end

% Set up global plotting settings for consistency
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultTextInterpreter', 'tex');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

% Create figure with specified dimensions
fig_handle = figure;
set(fig_handle, 'Units', 'centimeters');
set(fig_handle, 'Position', [0 0 fig_width fig_height]);

% Get coordinate conversion factors
lat_center = mean(terrain_data.y);
lon_scale = 111000 * cos(lat_center * pi/180);  % meters per degree longitude
lat_scale = 111000;  % meters per degree latitude

% Convert to meters
x_meters = (terrain_data.x - mean(terrain_data.x)) * lon_scale;
y_meters = (terrain_data.y - mean(terrain_data.y)) * lat_scale;

% Create terrain mesh
[X_terrain, Y_terrain] = meshgrid(x_meters, y_meters);

% Plot terrain surface
terrain_handle = surf(X_terrain, Y_terrain, terrain_data.terrain, 'EdgeColor', 'none', 'FaceAlpha', 1);

% Create continuous green to brown colormap (more gradual transition)
num_colors = 256;  % More colors for smoother gradient
green_low = [0.2, 0.6, 0.2];   % Dark green (low elevation)
green_mid = [0.5, 0.8, 0.3];   % Light green (mid-low elevation)
yellow = [0.8, 0.8, 0.3];      % Yellow (middle elevation)
brown_light = [0.7, 0.6, 0.3]; % Light brown (mid-high elevation)
brown_dark = [0.5, 0.3, 0.1];  % Dark brown (high elevation)

% Create interpolation points (0 to 1)
positions = [0, 0.25, 0.5, 0.75, 1];
colors = [green_low; green_mid; yellow; brown_light; brown_dark];

% Interpolate to create smooth gradient
r = interp1(positions, colors(:,1), linspace(0,1,num_colors), 'pchip');
g = interp1(positions, colors(:,2), linspace(0,1,num_colors), 'pchip');
b = interp1(positions, colors(:,3), linspace(0,1,num_colors), 'pchip');

terrain_colormap = [r' g' b'];

% Apply custom colormap
colormap(gca, terrain_colormap);

% Add colorbar
c = colorbar;
ylabel(c, 'Elevation (m)', 'FontWeight', 'bold');

% Enhance lighting
lighting gouraud
material dull
light('Position', [0.1 0.4 0.4], 'Style', 'infinite');
% light('Position', [-0.3, 0.2, 0.6], 'Style', 'infinite');

% Format axes and labels
xlabel('East (m)', 'FontWeight', 'bold');
ylabel('North (m)', 'FontWeight', 'bold');
zlabel('Elevation (m)', 'FontWeight', 'bold');

% Calculate plot area size from figure size and padding
plot_width = fig_width - left_padding - right_padding;
plot_height = fig_height - bottom_padding - top_padding;

% Set position of axes in the figure
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [left_padding bottom_padding plot_width plot_height]);

% Set view and grid
grid off;
view(view_azimuth, view_elevation);

% Set axis limits and scaling
axis tight;
axis equal;

% Box on for cleaner appearance
box on;

% Add title
% title('Wicklow Mountains Terrain', 'FontWeight', 'bold');

% Print elevation range for reference
elev_min = min(min(terrain_data.terrain));
elev_max = max(max(terrain_data.terrain));
fprintf('Elevation range: %.1f to %.1f meters\n', elev_min, elev_max);

% Print completion message
fprintf('Terrain visualization complete\n');