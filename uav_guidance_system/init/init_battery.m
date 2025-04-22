function battery = init_battery()
    % INIT_BATTERY Initialize battery parameters
    %
    % A simple battery model for tracking energy consumption based on thrust levels
    % Returns a structure with battery properties
    
    % % Battery capacity in Watt-hours
    % battery.capacity_Wh = 200;  % TODO: Adjust based on actual UAV specifications
    % 
    % % Basic battery consumption model parameters
    % battery.base_power = 10;     % Power for avionics/systems at idle (W)
    % battery.thrust_coefficient = 0.5;  % Power per unit thrust squared (W/N²)

    
    % INIT_BATTERY Initialise battery parameters for a delivery drone
    %
    % Returns a structure with battery properties for energy consumption tracking
    % Tailored for a 9kg fixed-wing delivery drone with 3m wingspan
    
    % Battery capacity calculation based on typical energy density
    % Using 180 Wh/kg for modern LiPo batteries with packaging overhead
    battery_weight_kg = 0.1; % ---> 1.8; % Assuming battery is ~20% of aircraft mass
    battery.capacity_Wh = battery_weight_kg * 180; % Gives approximately 324 Wh
    
    % Basic battery consumption model parameters
    battery.base_power = 1; % ----> 25; % Higher base power for delivery systems (W)
    
    % Power calculations based on typical efficiency metrics for delivery drones
    % For a 9kg aircraft, expect ~40-50W per kg at cruise
    % Using thrust-to-weight ratio of approximately 0.4 for cruise
    % This gives nominal thrust of 9kg*9.81*0.4 ≈ 35N
    % And cruise power of ~400W, so ~11.4W/N
    battery.thrust_coefficient = 12; % Power per unit thrust (W/N)
    
    % For more realistic simulation
    battery.discharge_efficiency = 0.85; % Battery discharge efficiency
    battery.reserve_capacity = 0.2; % 20% safety margin
    battery.usable_capacity_Wh = battery.capacity_Wh * (1 - battery.reserve_capacity);
    
    % For information display
    battery.initial_capacity_Wh = battery.capacity_Wh; % starting on low battery to see regen

end
