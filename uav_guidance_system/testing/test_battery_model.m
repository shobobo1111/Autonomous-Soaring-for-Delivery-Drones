% TEST_BATTERY Simple test of the battery model
%
% This script simulates a flight with different thrust levels
% and calculates the energy consumption

% Clear workspace
clear;
clc;

% Initialize battery
battery = init_battery();

% Simulate a flight with different thrust levels
% [0, 5, 15] = [gliding, low, full power]
thrust_history = [
    15, 15, 15, 15, 15,  % Takeoff & climb (full power)
    5, 5, 5, 5, 5,       % Cruise (low power)
    0, 0, 0, 0, 0,       % Glide (no power)
    5, 5, 5, 5, 5        % Approach (low power)
];

% Time steps (seconds)
dt = 10; % 10 seconds per step
time_history = (0:length(thrust_history)) * dt;  % One more time point than thrust points

% Calculate energy usage
energy_data = calculate_energy_usage(battery, thrust_history, time_history);

% Display results
disp('Flight Summary:');
disp(['Duration: ', num2str(time_history(end)), ' seconds']);
disp(['Energy consumed: ', num2str(energy_data.total_energy_Wh), ' Wh']);
disp(['Battery remaining: ', num2str(energy_data.percent_remaining), '%']);

% Plot thrust and power over time
figure;
subplot(2,1,1);
% Use stairs to show the thrust level held constant during each step
stairs(time_history(1:end-1), thrust_history, 'LineWidth', 2);
title('Thrust Level Over Time');
xlabel('Time (s)');
ylabel('Thrust');
grid on;

subplot(2,1,2);
% Plot power consumption over time
stairs(time_history(1:end-1), energy_data.power_values, 'LineWidth', 2);
title('Power Consumption Over Time');
xlabel('Time (s)');
ylabel('Power (W)');
grid on;

% TODO: Add more visualizations for energy consumption