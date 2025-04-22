% generate_multi_step_path.m
function complete_path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery, varargin)
    % Generates a path consisting of multiple sequential optimised steps
    %
    % Inputs:
    %   initial_state: Initial aircraft state
    %   aircraft: Aircraft parameters
    %   wind_field: Wind field structure
    %   num_steps: Number of optimization steps to take
    %
    % Output:
    %   complete_path: Structure containing full path information

        % Check for time_to_deliver parameter and pass to init_mission if provided
    time_to_deliver = [];
    if nargin > 5 && ~isempty(varargin{1})
        time_to_deliver = varargin{1};
        fprintf('Using custom time_to_deliver: %.1f seconds\n', time_to_deliver);
    end
    
    % Load mission parameters with custom time if provided
    mission = init_mission(time_to_deliver);
    
    
    % Initialise storage for complete path
    complete_path.states.x = [];
    complete_path.states.y = [];
    complete_path.states.z = [];
    complete_path.states.V = []; 
    complete_path.states.V_rel = []; 
    complete_path.states.phi = []; % Bank angle history
    complete_path.states.psi = []; % heading angle history
    complete_path.states.gamma = []; % pitch angle history
    complete_path.rewards = [];
    complete_path.states.thrust = [];
    % complete_path.steps = [];
    complete_path.states.mode = []; % Added: Operating mode history 12/04/2025
    complete_path.states.battery_per = []; % Added: Battery percentage history 12/04/2025

    elapsed_time = 0;
    current_battery = battery.capacity_Wh * 0.5; % was set at  * 0.5 for short tests
    
    % Start from initial state
    current_state = initial_state;

    % load in waypoint position for termination condition
    waypoint = init_waypoint();
    proximity_to_waypoint = sqrt((waypoint.x - current_state.x)^2 + (waypoint.y - current_state.y)^2 ); % updated this to consider arrival to waypoint excluding vertical component as a package delivery could be made from many heights

    arrival_tol = waypoint.arrival_tol;

    % loading in mission parameters
    % mission = init_mission();

    % Initialise step counter
    step = 1;

    %bool variable to stop run if in 
    reached_waypoint = false;
    
    % For each step
    % for step = 1:num_steps
    while step <= num_steps && current_battery > 0 && ~reached_waypoint
    % while step <= num_steps && current_battery > 0 % TEMPORARILY REMOVE PROX TO WAYPOINT ^ just for loiter opt

        % fprintf('Step %d starting at | x: %.0fm | y: %.0fm | z :%.0fm | at time %.1f seconds\n', ...
        %     step, current_state.x, current_state.y, current_state.z, elapsed_time);
        % fprintf('Angles: phi=%.1f°, gamma=%.1f°, psi=%.1f°\n\n', ...
        %         current_state.phi, current_state.gamma, current_state.psi);

        % Calculate current battery percentage
        battery_percentage = (current_battery / battery.initial_capacity_Wh) * 100;

        % Determine operating mode based on current conditions
        % =================== MODE DETERMINATION LOGIC ===================
        % Default mode is balanced (constrained time)
        current_mode = 2; % should be 5

        % if proximity_to_waypoint < waypoint.loiter_within_r
        %     current_mode = 1; % loitering
        % end
        


        % % Time-based considerations
        % if elapsed_time < 10
        %     current_mode = 3; % Long distance mode at the start
        % end
        % 
        % % % Battery-based considerations
        % % if battery_percentage < 30
        % %     current_mode = 5; % Energy saving mode when battery is low
        % % end
        % % 
        % % Altitude-based considerations
        % if current_state.z < mission.alt_floor + mission.alt_floor_tol
        %     current_mode = 2; % Urgent/straight to point when altitude is low
        % end

        % wind_diff = abs(current_state.V - current_state.V_rel);
        % fprintf('V - V_rel Difference :%.5f\n', wind_diff);
        % 
        % %% IT IS UNCERTAIN IF THIS WORKS WITH RIDGE LIFT
        % if wind_diff > 0.01
        %     current_mode = 4; 
        % end


        wind = get_wind_at_position(current_state.x, current_state.y, current_state.z, wind_field);
        fprintf('%2.2f',wind.U);
    

        % if wind.U > 1
        %     current_mode = 4; 
        % end


        
        % % Proximity-based considerations
        % if proximity_to_waypoint < 100
        %     current_mode = 2; % Urgent/straight to point when close to waypoint 
        % end
        
        % Lateness considerations
        given_time = mission.time_to_deliver;
        time_left = given_time - elapsed_time;
        predicted_time = proximity_to_waypoint / aircraft.cruise_speed;
        p_time_w_buffer = predicted_time + (predicted_time * 0.5);

        if proximity_to_waypoint > waypoint.loiter_within_r
            lateness = p_time_w_buffer - time_left;
        
            if lateness > 0
                current_mode = 2; % Urgent mode if we're going to be late
            end
        else
            lateness = -65;
        end


        % 
        % % Waypoint handling
        % if proximity_to_waypoint < waypoint.loiter_within_r
        %     current_mode = 1; % Loiter mode if near waypoint but still high
        % end
        % 
        % Special battery regeneration case
        % Conditions: good altitude, not urgent timing, decent battery left
        % if current_state.z > (waypoint.z+50) && lateness < -60 && battery_percentage > 30 && battery_percentage < 99.9
        %     current_mode = 6; % Regeneration mode
        % end
        % =============== END MODE DETERMINATION LOGIC =================


            % Determine if waypoint is reached based on multiple conditions
        if waypoint.loiter_within_r > 0
            % If loiter_within_r is set, aircraft must stay within this radius
            % (This will keep reached_waypoint as false, forcing the aircraft to loiter)
            reached_waypoint = false;
        else
            % Otherwise, waypoint is reached if:
            % 1. Aircraft is within arrival tolerance horizontally
            % 2. Aircraft is above the waypoint altitude
            reached_waypoint = (proximity_to_waypoint <= waypoint.arrival_tol) && ...
                              (current_state.z >= (waypoint.z - waypoint.arrival_tol));
        end


        % Improved formatted output for each step
        fprintf('\n====== Step %d (%0.1f sec) ======\n', step, elapsed_time);
        fprintf('Position: [%5.0f, %5.0f, %5.0f] m  V_rel: %5.4fm/s V: %5.4fm/s\n', current_state.x, current_state.y, current_state.z, current_state.V_rel, current_state.V);


        % Display current mode
        mode_names = {'Loiter', 'Urgent', 'Long Distance', 'Updraft Climb', 'Energy Saving', 'Regeneration'};
        fprintf('Operating mode: %s\n', mode_names{current_mode});

        % Generate and evaluate possible paths from current position
        paths = generate_paths(current_state, aircraft, wind_field, elapsed_time, current_mode);

        % plot_path_flower(paths, wind_field, step, elapsed_time);  % plots the flower at each point
        
        % Best path is first one (already sorted by reward in generate_paths)
        best_path = paths(1);
        
        % % Store this segment
        % complete_path.states.x = [complete_path.states.x, best_path.states.x];
        % complete_path.states.y = [complete_path.states.y, best_path.states.y];
        % complete_path.states.z = [complete_path.states.z, best_path.states.z];
        % complete_path.states.V = [complete_path.states.V, best_path.states.V];
        % complete_path.states.V_rel = [complete_path.states.V_rel, best_path.states.V_rel];

        % Modified code to skip the first point of subsequent segments (except for the first segment)
        if step == 1
            % For the first step, include the entire path
            complete_path.states.x = [complete_path.states.x, best_path.states.x];
            complete_path.states.y = [complete_path.states.y, best_path.states.y];
            complete_path.states.z = [complete_path.states.z, best_path.states.z];
            complete_path.states.V = [complete_path.states.V, best_path.states.V];
            complete_path.states.V_rel = [complete_path.states.V_rel, best_path.states.V_rel];
        else
            % For subsequent steps, skip the first point to avoid duplication
            complete_path.states.x = [complete_path.states.x, best_path.states.x(2:end)];
            complete_path.states.y = [complete_path.states.y, best_path.states.y(2:end)];
            complete_path.states.z = [complete_path.states.z, best_path.states.z(2:end)];
            complete_path.states.V = [complete_path.states.V, best_path.states.V(2:end)];
            complete_path.states.V_rel = [complete_path.states.V_rel, best_path.states.V_rel(2:end)];
        end

        complete_path.states.phi = [complete_path.states.phi, best_path.bank_angle]; 
        complete_path.states.gamma = [complete_path.states.gamma, best_path.pitch_angle];
        complete_path.states.psi = [complete_path.states.psi, current_state.psi]; % need to get this from flight_dynamics calcs
        complete_path.states.thrust = [complete_path.states.thrust, best_path.thrust];
        complete_path.rewards = [complete_path.rewards, best_path.total_reward];
        complete_path.states.mode = [complete_path.states.mode, current_mode]; % Store mode
        complete_path.states.battery_per = [complete_path.states.battery_per, battery_percentage]; % Store battery 
        
        % checking how close UAV is to waypoint to inform while loop exit
        proximity_to_waypoint = sqrt((waypoint.x - current_state.x)^2 + ...
                (waypoint.y - current_state.y)^2); % updated this to consider arrival to waypoint excluding vertical component as a package delivery could be made from many heights
            %+ (waypoint.z - current_state.z)^2
                
        % % Basic battery consumption model parameters
        % base_power = battery.base_power;     % Power for avionics/systems at idle (W)
        % thrust_coefficient = battery.thrust_coefficient;  % Power per unit thrust squared (W/N²)
        % 
        % thrust_level = best_path.thrust;
        % power_consumed_this_step = base_power + thrust_coefficient*thrust_level^2; 
        % dt = 2; % 2 second time step, remember to change this if changed else where
        % energy_consumed_this_step = power_consumed_this_step * dt;
        % current_battery = current_battery - energy_consumed_this_step;
        % fprintf('Thrust level: %.2f \nCurrent Battery Remaining: %.2f \nBattery used: %.2f\n',thrust_level,current_battery,energy_consumed_this_step);        
        % 
        %% Battery tracking
        % Calculate energy consumption for this step
        base_power = battery.base_power; % Power for systems (W)
        thrust_coefficient = battery.thrust_coefficient; % Power per unit thrust (W/N)
        thrust_level = best_path.thrust; % Current thrust in N
        
        % Calculate power consumption (W)
        power_consumed_this_step = base_power + thrust_coefficient * thrust_level;
        
        % Account for inefficiency in power delivery
        power_consumed_this_step = power_consumed_this_step / battery.discharge_efficiency;
        
        % Calculate energy consumed in this step (Wh)
        dt = 2; % 2 second time step
        energy_consumed_this_step = power_consumed_this_step * (dt/3600); % Convert to Wh
        
        % Update remaining battery capacity
        current_battery = current_battery - energy_consumed_this_step;
        
        % % Calculate percentage remaining for display
        % percent_remaining = (current_battery / battery.initial_capacity_Wh) * 100;
        % 
        % fprintf('Thrust: %.1f N | Power: %.1f W | Battery: %.1f Wh (%.1f%%) | Used this step: %.3f Wh\n\n', thrust_level, power_consumed_this_step, current_battery, percent_remaining, energy_consumed_this_step);


        % Battery status with visual indicator
        percent_remaining = (current_battery / battery.initial_capacity_Wh) * 100;
        battery_bar = repmat('█', 1, round(percent_remaining/5)); % Each block = 5%
        battery_empty = repmat('░', 1, 20-round(percent_remaining/5));
        fprintf('Battery: %s%s %.1f%% (%.1f Wh)\n', battery_bar, battery_empty, percent_remaining, current_battery);

        % Thrust and power metrics
        fprintf('Thrust: %3.1f N | Power: %5.1f W | Used: %.3f Wh\n', thrust_level, power_consumed_this_step, energy_consumed_this_step);


        %% Lateness, Mission time and waypoint progress

        % Calculate navigation reward with time pressure
        % d_current = sqrt((waypoint.x - current_state.x)^2 + ...
        %                 (waypoint.y - current_state.y)^2 + ...
        %                 (waypoint.z - current_state.z)^2);
        % d_next = sqrt((waypoint.x - next_state.x)^2 + ...
        %               (waypoint.y - next_state.y)^2 + ...
        %               (waypoint.z - next_state.z)^2);
                      
        % given_time = mission.time_to_deliver; % the time allocated for the delivery
        % time_left = given_time - elapsed_time; % given 100s - 20s past already -> you have 80s left
        % 
        % % the predicted time it will take to get to the waypoint at cruise
        % % added buffer to warn of lateness before it happens
        % predicted_time = proximity_to_waypoint / aircraft.cruise_speed; % 150m left to go / 15m/s -> 10s predicted time to get there
        % p_time_w_buffer = predicted_time + (predicted_time * 0.5); % now 15s to get there cause life
        % 
        % % lateness is how late you would arrive at cruise (plus buffer)
        % lateness = p_time_w_buffer - time_left; % if predicted time is 35 and time left is 30, lateness is 5

        % fprintf('Given Time: %.1f s | Time Left: %.1f s | lateness: %.1f  | \n\n', given_time, time_left,lateness);

        % Waypoint navigation status
        fprintf('Distance to waypoint: %5.1f m\n', proximity_to_waypoint);
        
        % Mission timing status
        remaining_bar_width = 20;
        time_progress = min(100, (elapsed_time / given_time) * 100);
        time_bar = repmat('█', 1, round(time_progress/5));
        time_remaining = repmat('░', 1, remaining_bar_width-round(time_progress/5));
        fprintf('Mission time: %s%s %.1f%% (%0.1f/%0.1f sec)\n', time_bar, time_remaining, time_progress, elapsed_time, given_time);
        
        % Add warning indicator for lateness
        if lateness > 0
            fprintf('⚠️  LATENESS WARNING: Projected to arrive %.1f seconds late\n', lateness);
        end

            
        
        
        %% Updating Positions & angles
        % Update current state for next step using final state of best path
        % Get final index
        final_idx = length(best_path.states.x);
        
        % Update position
        current_state.x = best_path.states.x(final_idx);
        current_state.y = best_path.states.y(final_idx);
        current_state.z = best_path.states.z(final_idx);
        current_state.V = best_path.states.V(final_idx);
        current_state.V_rel = best_path.states.V_rel(final_idx);
        
        % Update angles
        current_state.phi = best_path.bank_angle;     % Bank angle from chosen path
        current_state.gamma = best_path.pitch_angle;  % Pitch angle from chosen path
        current_state.psi = best_path.states.psi(final_idx);
        
        % % Calculate new heading based on path taken
        % if final_idx > 1
        %     % dx = best_path.states.x(final_idx) - best_path.states.x(final_idx-1);
        %     % dy = best_path.states.y(final_idx) - best_path.states.y(final_idx-1);
        %     current_state.psi = atan2(dy, dx);  % New heading based on actual motion
        % end

        % plot_path_flower(paths, wind_field, step, elapsed_time)
        elapsed_time = elapsed_time + 2; % 2 is the predicition time
        step = step + 1;
    end
    
    complete_path.steps = step - 1;
    % Determine which condition triggered the termination
    % Determine which condition triggered the termination
    if current_battery <= 0
        termination_reason = "Battery depleted :/";
    elseif reached_waypoint
        termination_reason = "Waypoint reached!";
    else
        termination_reason = "Maximum steps reached";
    end

    fprintf('\n___________________________\nMission ended: %s at step %d\n\n', termination_reason, complete_path.steps);
end

