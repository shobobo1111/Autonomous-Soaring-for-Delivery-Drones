function energy_data = calculate_energy_usage(battery, thrust_history, time_history)
    % CALCULATE_ENERGY_USAGE Calculate the energy usage based on thrust history
    %
    % Inputs:
    %   battery: Battery structure from init_battery()
    %   thrust_history: Array of thrust values used during flight
    %   time_history: Array of timestamps corresponding to thrust values
    %
    % Outputs:
    %   energy_data: Structure with energy consumption data
    
    % Check input sizes
    if length(thrust_history) + 1 ~= length(time_history)
        error('Time history should have one more entry than thrust history');
    end
    
    % Get the time differences between steps
    if length(time_history) <= 1
        energy_data.total_energy_Wh = 0;
        energy_data.percent_remaining = 100;
        return;
    end
    
    % Map thrust values to power consumption
    thrust_levels = [0, 5, 15]; % Same as in generate_paths
    power_values = zeros(size(thrust_history));
    
    for i = 1:length(thrust_history)
        % Find which thrust level was used
        [~, idx] = min(abs(thrust_levels - thrust_history(i)));
        power_values(i) = battery.power_consumption(idx);
    end
    
    % Calculate time differences for each step
    dt = diff(time_history);  % This gives us the duration of each step
    
    % Calculate energy usage for each time step (Wh)
    % Power (W) * time (h) = energy (Wh)
    energy_per_step = zeros(size(thrust_history));
    for i = 1:length(thrust_history)
        energy_per_step(i) = power_values(i) * (dt(i) / 3600);  % Convert seconds to hours
    end
    
    % Calculate total energy consumed
    total_energy_Wh = sum(energy_per_step);
    
    % Update battery state
    remaining_charge = battery.capacity_Wh - total_energy_Wh;
    percent_remaining = (remaining_charge / battery.capacity_Wh) * 100;
    
    % Ensure we don't go below 0%
    if percent_remaining < 0
        percent_remaining = 0;
        remaining_charge = 0;
    end
    
    % Prepare return data
    energy_data.total_energy_Wh = total_energy_Wh;
    energy_data.power_values = power_values;
    energy_data.energy_per_step = energy_per_step;
    energy_data.percent_remaining = percent_remaining;
    energy_data.remaining_charge_Wh = remaining_charge;
end