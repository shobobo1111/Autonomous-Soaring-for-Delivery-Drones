terrain_data = load('terrain_data.mat');

full_path = fullfile('simulation_data', '20250421', '20_flight_data_12:26.mat');

% plot_and_save('', '', '', '', '', full_path); 
plot_and_save('', '', '', '', terrain_data, full_path);