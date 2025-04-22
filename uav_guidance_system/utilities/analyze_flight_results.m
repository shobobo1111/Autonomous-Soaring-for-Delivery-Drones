% analyze_flight.m - Script to analyze loaded flight data

% Check if flight data is loaded
if ~exist('complete_path', 'var')
    error('No flight data found. Please load a .mat file first.')
end

% Load battery parameters
battery = init_battery();

% BATTERY USAGE
thrust_levels = complete_path.states.thrust;
base_power = battery.base_power; 
thrust_coefficient = battery.thrust_coefficient;
dt = 2; % seconds per step

% Calculate energy consumption for each step - using exact method from your flight code
total_energy_consumed = 0;
for i = 1:length(thrust_levels)
    % Power in Watts
    power_consumed = base_power + thrust_coefficient * thrust_levels(i);
    
    % Account for inefficiency in power delivery
    power_consumed = power_consumed / battery.discharge_efficiency;
    
    % Energy in Watt-hours
    energy_consumed = power_consumed * (dt/3600);
    total_energy_consumed = total_energy_consumed + energy_consumed;
end

% Calculate battery percentages
initial_battery = battery.capacity_Wh;
final_battery = initial_battery - total_energy_consumed;
battery_used_percent = (total_energy_consumed / initial_battery) * 100;

% ALTITUDE METRICS
initial_altitude = initial_state.z;
max_altitude = max(complete_path.states.z);
final_altitude = complete_path.states.z(end);
altitude_gain = max_altitude - initial_altitude;

% TIME METRICS
flight_time = complete_path.steps * 2; % seconds

% WAYPOINT DISTANCE
waypoint = init_waypoint();
final_position = [complete_path.states.x(end), complete_path.states.y(end), complete_path.states.z(end)];
waypoint_position = [waypoint.x, waypoint.y, waypoint.z];
final_distance = norm(final_position - waypoint_position);

% MISSION TIMING
mission = init_mission();
time_to_deliver = mission.time_to_deliver;

% Calculate delivery status
if flight_time > time_to_deliver
    delivery_status = sprintf('Late by %.1f seconds', flight_time - time_to_deliver);
elseif final_distance > waypoint.arrival_tol
    delivery_status = sprintf('Incomplete (%.1fm from waypoint)', final_distance);
elseif flight_time < time_to_deliver
    delivery_status = sprintf('Early by %.1f seconds', time_to_deliver - flight_time);
else
    delivery_status = 'On time';
end

% OPERATING MODE STATISTICS
if isfield(complete_path.states, 'mode')
    modes = complete_path.states.mode;
    mode_names = {'Loiter', 'Urgent', 'Long Distance', 'Balanced', 'Energy Saving', 'Regen'};
    
    % Calculate time in each mode
    mode_times = zeros(1, length(mode_names));
    for i = 1:length(mode_names)
        mode_times(i) = sum(modes == i) * 2; % seconds
    end
    
    % Calculate percentages
    mode_percentages = (mode_times / flight_time) * 100;
else
    % If modes not available
    mode_names = {'Loiter', 'Urgent', 'Long Distance', 'Balanced', 'Energy Saving', 'Regen'};
    mode_times = zeros(1, 6);
    mode_times(2) = flight_time; % Assume all urgent by default
    mode_percentages = (mode_times / flight_time) * 100;
end

% PRINT RESULTS
fprintf('\n----- FLIGHT ANALYSIS RESULTS -----\n');
fprintf('Battery usage: %.1f Wh (%.1f%%)\n', total_energy_consumed, battery_used_percent);
fprintf('Initial capacity: %.1f Wh\n', initial_battery);
fprintf('Remaining capacity: %.1f Wh\n\n', final_battery);

fprintf('Max altitude gain: %.1f meters\n', altitude_gain);
fprintf('Initial altitude: %.1f meters\n', initial_altitude);
fprintf('Final altitude: %.1f meters\n\n', final_altitude);

fprintf('Total flight time: %.1f seconds\n', flight_time);
fprintf('Delivery status: %s\n', delivery_status);
fprintf('Distance to waypoint: %.1f meters\n\n', final_distance);

fprintf('Time spent in each mode:\n');
for i = 1:length(mode_names)
    if mode_times(i) > 0
        fprintf('  %s: %.1f seconds (%.1f%%)\n', mode_names{i}, mode_times(i), mode_percentages(i));
    end
end