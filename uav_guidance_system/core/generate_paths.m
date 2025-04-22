function [paths] = generate_paths(current_state, aircraft, wind_field, elapsed_time, operating_mode)
    % GENERATE_PATHS Creates and evaluates possible flight paths from current state
    % including "feeler lines" that extend beyond each path to evaluate future potential
    
    % Get configs
    feeler_config = init_feeler_config();
    waypoint = init_waypoint();
    
    % Define control options for path generation
    bank_angles = linspace(-60, 60, 7) * (pi/180);   % (bottom, top, number) -> Converted to radians 
    % pitch_angles = linspace(-1.43239, 15, 4) * (pi/180);  
    pitch_angles = [-1.43239, 5, 15] * (pi/180);  % In radians
    thrust_levels = [0, 5, 10];  % [gliding, low power, full power] (minus representing regenerative braking)

    wind_diff = abs(current_state.V - current_state.V_rel);
    fprintf('V - V_rel Difference :%.5f\n', wind_diff);

    wind = get_wind_at_position(current_state.x, current_state.y, current_state.z, wind_field);
    if operating_mode == 4
        pitch_angles = [(-3*wind.U), -1.43239, 5, 15] * (pi/180);
    end

    % wind = get_wind_at_position(current_state.x, current_state.y, current_state.z, wind_field);
    % fprintf('%2.2f',wind.U);
    % 
    % if operating_mode ~= 2 && operating_mode ~= 6
    %     if wind.U > 2
    %         pitch_angles = [-10, -5, -1.43239, 5, 15] * (pi/180); % more angles for higher updraft speeds
    %         operating_mode = 4; % set operating mode for plotting
    %     elseif wind.U > 1
    %         pitch_angles = [-5, -1.43239, 5, 15] * (pi/180); 
    %         operating_mode = 4;
    %     else
    %         pitch_angles = [-1.43239, 5, 15] * (pi/180); 
    %     end
    % end

    % if operating_mode == 4 
    %     pitch_angles = [-10, -5, -1.43239, 5, 15] * (pi/180);  % In radians
    % else 
    %     pitch_angles = [-1.43239, 5, 15] * (pi/180);  % In radians
    % end
    
    

    if operating_mode == 6
        thrust_levels = [-2, -1, 0]; 
        pitch_angles = [(-3*wind.U), -1.43239, 5, 15] * (pi/180);
    else
        thrust_levels = [0, 5, 10]; 
    end
    
    % Calculate total number of possible paths
    num_paths = length(bank_angles) * length(pitch_angles) * length(thrust_levels);
    
    % Initialise paths structure array
    paths = struct('states', cell(1, num_paths), ...
                  'feeler', cell(1, num_paths), ...    
                  'rewards', cell(1, num_paths));
    
    % Time settings for path prediction
    prediction_time = 2;    % Length of main path prediction (seconds)
    num_steps = 20;        % Number of discrete steps in prediction
    dt = prediction_time/num_steps;
    
    % Generate path for each control combination
    path_idx = 1;
    for i_bank = 1:length(bank_angles)
        for i_pitch = 1:length(pitch_angles)
            for i_thrust = 1:length(thrust_levels)

                % initilaising a flag to track if paths are safe to take
                % paths(path_idx).unsafe_path = false;

                % Store control parameters for this path
                paths(path_idx).bank_angle = bank_angles(i_bank);
                paths(path_idx).pitch_angle = pitch_angles(i_pitch);
                paths(path_idx).thrust = thrust_levels(i_thrust);
                
                % this initialises the state storage arrays
                paths(path_idx).states.x = zeros(1, num_steps+1);
                paths(path_idx).states.y = zeros(1, num_steps+1);
                paths(path_idx).states.z = zeros(1, num_steps+1);
                paths(path_idx).states.V = zeros(1, num_steps+1);
                paths(path_idx).states.V_rel = zeros(1, num_steps+1);
                
                
                % Set control inputs for dynamics
                control_input.phi_command = bank_angles(i_bank);
                control_input.gamma_command = pitch_angles(i_pitch);
                control_input.thrust = thrust_levels(i_thrust);
                
                % Initialise reward storage
                paths(path_idx).rewards = zeros(1, num_steps);
                
                % Start from current state
                current = current_state;
                
                % Store initial conditions
                paths(path_idx).states.x(1) = current.x;
                paths(path_idx).states.y(1) = current.y;
                paths(path_idx).states.z(1) = current.z;
                paths(path_idx).states.V(1) = current.V;
                paths(path_idx).states.V_rel(1) = current.V_rel;
                
                % Calculate main path trajectory and rewards
                total_reward = 0;
                for j = 1:num_steps
                    local_wind = get_wind_at_position(current.x, current.y, current.z, wind_field);
                    previous = current; 
                    current = flight_dynamics(current, control_input, dt, local_wind, aircraft);

                    %     % Immediate safety check - exit if angle exceeds limit or speed is too low
                    % if abs(current.gamma) > aircraft.max_climb_angle % || current.V_rel < (aircraft.min_speed + 1.0) 
                    %     paths(path_idx).unsafe_path = true;
                    %     fprintf('a path invalidated')
                    % end
                    
                    % Calculate and store reward
                    step_reward = calculate_reward(previous, current, local_wind, waypoint, elapsed_time, control_input, operating_mode);
                    paths(path_idx).rewards(j) = step_reward;
                    total_reward = total_reward + step_reward;
                    
                    % Store states
                    paths(path_idx).states.x(j+1) = current.x;
                    paths(path_idx).states.y(j+1) = current.y;
                    paths(path_idx).states.z(j+1) = current.z;
                    paths(path_idx).states.V(j+1) = current.V;
                    paths(path_idx).states.V_rel(j+1) = current.V_rel;
                    paths(path_idx).states.psi(j+1) = current.psi;  % Add this line to store heading
                end
                
                if feeler_config.enabled
                    % Calculate final heading and pitch for feeler line
                    final_idx = length(paths(path_idx).states.x);
                    if final_idx > 1
                        % Calculate heading 
                        final_heading = atan2(paths(path_idx).states.y(end) - paths(path_idx).states.y(end-1), ...
                                           paths(path_idx).states.x(end) - paths(path_idx).states.x(end-1));
                
                        % Calculate pitch (vertical angle)
                        final_pitch = atan2(paths(path_idx).states.z(end) - paths(path_idx).states.z(end-1), ...
                        sqrt((paths(path_idx).states.x(end) - paths(path_idx).states.x(end-1))^2 + ...
                        (paths(path_idx).states.y(end) - paths(path_idx).states.y(end-1))^2));
                    else
                        final_heading = current_state.psi;
                        final_pitch = current_state.gamma;
                    end
                
                    aircraft.max_climb_angle;  % 15 degrees max climb angle
    
                    if abs(final_pitch) <= aircraft.max_climb_angle % 15 degrees max climb angle
                        %     paths(path_idx).feeler_reward = 0;
                        %     continue;  % Skip feeler line generation for this path
                        % end
                    
                        % Generate feeler line points maintaining heading and pitch
                        dx = feeler_config.length * cos(final_pitch) * cos(final_heading);
                        dy = feeler_config.length * cos(final_pitch) * sin(final_heading); 
                        dz = feeler_config.length * sin(final_pitch);
                        
                        % Generate feeler line points
                        paths(path_idx).feeler.x = linspace(paths(path_idx).states.x(end), ...
                            paths(path_idx).states.x(end) + dx, ...
                            feeler_config.num_points);
                        paths(path_idx).feeler.y = linspace(paths(path_idx).states.y(end), ...
                            paths(path_idx).states.y(end) + dy, ...
                            feeler_config.num_points);
                        paths(path_idx).feeler.z = linspace(paths(path_idx).states.z(end), ...
                            paths(path_idx).states.z(end) + dz, ...
                            feeler_config.num_points);
                        
                        % Calculate feeler rewards
                        feeler_reward = 0;
                        for f = 1:feeler_config.num_points-1
                            % Create state structures for reward calculation
                            feeler_state = struct('x', paths(path_idx).feeler.x(f), ...
                             'y', paths(path_idx).feeler.y(f), ...
                             'z', paths(path_idx).feeler.z(f), ...
                             'V', paths(path_idx).states.V(end), ... % Use final velocity from main path
                             'V_rel', paths(path_idx).states.V_rel(end), ...
                             'psi', final_heading, ...               % Use calculated final heading
                             'gamma', final_pitch, ...               % Use calculated final pitch
                             'phi', paths(path_idx).bank_angle);     % Use path bank angle
        
                            next_feeler_state = struct('x', paths(path_idx).feeler.x(f+1), ...
                                  'y', paths(path_idx).feeler.y(f+1), ...
                                  'z', paths(path_idx).feeler.z(f+1), ...
                                  'V', paths(path_idx).states.V(end), ... % Maintain same velocity
                                  'V_rel', paths(path_idx).states.V_rel(end), ... % added to track V_rel for stall condition in calculate rewards
                                  'psi', final_heading, ...
                                  'gamma', final_pitch, ...
                                  'phi', paths(path_idx).bank_angle);
                            
                            % Get wind at feeler point
                            local_wind = get_wind_at_position(feeler_state.x, feeler_state.y, feeler_state.z, wind_field);
                            
                            % Calculate reward with decay and weight
                            step_reward = calculate_reward(feeler_state, next_feeler_state, local_wind, waypoint, elapsed_time, control_input, operating_mode);
                            feeler_reward = feeler_reward + step_reward * ...
                                feeler_config.initial_weight * ...
                                (feeler_config.decay_rate^f);
                        end
                        % Store feeler reward
                        paths(path_idx).feeler_reward = feeler_reward;
                        total_reward = total_reward + feeler_reward;
    
                    else
                        paths(path_idx).feeler_reward = 0;
                    end
                end
                % Store total reward
                paths(path_idx).total_reward = total_reward;
                path_idx = path_idx + 1;
                

                % %% checks if path is unsafe and assigns it -inf
                % if paths(path_idx).unsafe_path
                %     paths(path_idx).total_reward = -Inf;
                % else
                %     paths(path_idx).total_reward = total_reward;
                % end 
                % path_idx = path_idx + 1;
            end
        end
    end
    
    % Sort paths by total reward (best paths first)
    [~, reward_order] = sort([paths.total_reward], 'descend');
    paths = paths(reward_order);

end
