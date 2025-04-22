% test_wind_field.m
clear all; close all; clc;

% Define grid to match flight space
x_limits = [-100, 300];   % Plenty of room for paths to develop
y_limits = [-100, 300];   % Symmetric around start point
z_limits = [0, 200];      % From ground to above flight level
grid_spacing = 20;        % Coarse enough to see vectors clearly

% Create wind field structure
wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing);

% Add constant wind with thermal
thermal_center = [50, 50];
thermal_strength = 2;    % m/s vertical
thermal_radius = 20;     % meters
constant_wind = [0, 0.1, 0];  % [Wx, Wy, Wz] background wind (m/s)

% Fill wind components
for i = 1:length(wind_field.x)
    for j = 1:length(wind_field.y)
        for k = 1:length(wind_field.z)
            % Calculate distance from thermal center
            r = sqrt((wind_field.X(j,i,k) - thermal_center(1))^2 + ...
                    (wind_field.Y(j,i,k) - thermal_center(2))^2);
            
            % Combine thermal with constant wind
            wind_field.Wx(j,i,k) = constant_wind(1);
            wind_field.Wy(j,i,k) = constant_wind(2);
            wind_field.U(j,i,k) = thermal_strength * exp(-(r/thermal_radius)^2) + constant_wind(3);
        end
    end
end

% Plot full 3D wind field
plot_wind_field(wind_field);
