% optimize_reward_weights.m
clear all; close all; clc;

%% Define global variables for weights that calculate_reward.m will check
global REWARD_WEIGHTS ga_history date_folder;
REWARD_WEIGHTS = []; % Initially empty



%% Setup data collection for optimisation history
ga_history = struct();
ga_history.best_fitness = [];
ga_history.best_weights = [];
ga_history.generation = [];
ga_history.timestamp = datetime('now');

% Create base folder for data storage
base_folder = 'simulation_data';
if ~exist(base_folder, 'dir')
    mkdir(base_folder);
end

% Create date folder
current_date = datestr(now, 'yyyymmdd');
date_folder = fullfile(base_folder, current_date);
if ~exist(date_folder, 'dir')
    mkdir(date_folder);
end

% Create a unique filename for this optimization run
timestamp = datestr(now, 'HHMMSS');
optimization_id = sprintf('GA_opt_%s_%s', current_date, timestamp);




%% Define the optimisation problem
% We want to optimize the reward function weights:
% w_P (potential energy) - keep fixed
% w_K (kinetic energy) - keep fixed
% w_U (updraft utilization) - optimise
% w_N (navigation) - optimise
% w_S (stall prevention) - keep fixed
% w_T (thrust/energy penalty) - optimise
% w_F (floor breaching penalty) - keep fixed


% Define which weights to optimise (1 = optimise, 0 = keep fixed)
optimize_mask = [1, 1, 1, 1, 0, 1, 0]; % [w_P, w_K, w_U, w_N, w_S, w_T, w_F]

% % Define bounds for weights (lower and upper)
% lb = [1.2, 1.5, 0.5, 0.1, 3.0, 0.1];  % Lower bounds - fixed values for non-optimised weights
% ub = [1.2, 1.5, 5.0, 2.0, 3.0, 2.0];  % Upper bounds - fixed values for non-optimised weights

% Define bounds for weights (lower and upper)
lb = [0.1, 0.3, 1, 0.01, 5.0, 0.01, 1];  % Lower bounds - with 0.01 minimum
ub = [1.3, 1.9, 5.0, 2.0, 5.0, 1.0, 1];     % Upper bounds - unchanged

% Filter out only the weights we want to optimize
lb_opt = lb(optimize_mask == 1);
ub_opt = ub(optimize_mask == 1);
nvars = sum(optimize_mask);  % Number of variables to optimize

% Setup global variables to store our wind field and aircraft
global wind_field aircraft initial_state evaluation_waypoint original_weights optimize_indices battery;

% Store indices of weights to optimize for later reference
optimize_indices = find(optimize_mask == 1);

% Original weights for comparison
original_weights.w_P = 1.2;
original_weights.w_K = 1.5;
original_weights.w_U = 3.0;
original_weights.w_N = 0.01;
original_weights.w_S = 5.0;
original_weights.w_T = 0.1;
original_weights.w_F = 1;

% Load aircraft configuration
aircraft = init_aircraft();

mission = init_mission();

battery = init_battery();

%% Create wind field
% Define grid to match flight space
x_limits = [-700, 700];   
y_limits = [-700, 700];   
z_limits = [0, 800];     
grid_spacing = 20;   % 20 works for 500X500X500 space

% Create wind field structure
wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing);

% Add desired thermal pattern
wind_field = create_thermals(wind_field, 'simple_sys');

%% Initial conditions
initial_state.x = -500; % -100 for ridge
initial_state.y = -300; % 500 for ridge
initial_state.z = 200; % 700 for ridge
initial_state.V = 15;        % m/s
initial_state.V_rel = 15;    % m/s
initial_state.gamma = 0;     % rad
initial_state.phi = 0;       % rad
initial_state.psi = 0;       % rad

%% Get waypoint for evaluation
evaluation_waypoint = init_waypoint();

% %% Configure genetic algorithm options
% options = optimoptions('ga', ...
%     'PopulationSize', 20, ...         % Size of population - reduced for faster execution
%     'MaxGenerations', 12, ...         % Maximum number of generations
%     'FunctionTolerance', 1e-3, ...    % Tolerance for stopping
%     'PlotFcn', @gaplotbestf, ...      % Plot best fitness
%     'Display', 'iter', ...            % Show iteration information
%     'UseParallel', true, ...          % Use parallel computing if available
%     'MaxStallGenerations', 4);        % Stop if no improvement for 4 generations

%% Configure genetic algorithm options
options = optimoptions('ga', ...
    'PopulationSize', 20, ...         % Size of population - reduced for faster execution
    'MaxGenerations', 12, ...         % Maximum number of generations
    'FunctionTolerance', 1e-3, ...    % Tolerance for stopping
    'PlotFcn', @gaplotbestf, ...      % Plot best fitness
    'Display', 'iter', ...            % Show iteration information
    'UseParallel', true, ...          % Use parallel computing if available
    'MaxStallGenerations', 4, ...     % Stop if no improvement for 4 generations
    'OutputFcn', @saveGAData);        % Add our data collection function

%% Run the genetic algorithm
[optimal_weights_subset, final_fitness, exitflag, output] = ga(@evaluateWeights, nvars, [], [], [], [], lb_opt, ub_opt, [], options);

% Reconstruct full weights vector
optimal_weights_array = [original_weights.w_P, original_weights.w_K, original_weights.w_U, ...
                         original_weights.w_N, original_weights.w_S, original_weights.w_T, original_weights.w_F];
                     
% Update only the optimized weights
optimal_weights_array(optimize_indices) = optimal_weights_subset;

%% Display results
fprintf('\nOptimization Results:\n');
fprintf('Optimal weights found:\n');
fprintf('  w_P (potential energy): %.4f %s\n', optimal_weights_array(1), iif(optimize_mask(1), '', '(fixed)'));
fprintf('  w_K (kinetic energy): %.4f %s\n', optimal_weights_array(2), iif(optimize_mask(2), '', '(fixed)'));
fprintf('  w_U (updraft): %.4f %s\n', optimal_weights_array(3), iif(optimize_mask(3), '', '(fixed)'));
fprintf('  w_N (navigation): %.4f %s\n', optimal_weights_array(4), iif(optimize_mask(4), '', '(fixed)'));
fprintf('  w_S (stall prevention): %.4f %s\n', optimal_weights_array(5), iif(optimize_mask(5), '', '(fixed)'));
fprintf('  w_T (energy/thrust): %.4f %s\n', optimal_weights_array(6), iif(optimize_mask(6), '', '(fixed)'));
fprintf('  w_F (alt floor breach): %.4f %s\n', optimal_weights_array(7), iif(optimize_mask(7), '', '(fixed)'));
fprintf('Final fitness value: %.4f\n', final_fitness);

% Update optimal_weights structure for use in visualization
optimal_weights.w_P = optimal_weights_array(1);
optimal_weights.w_K = optimal_weights_array(2);
optimal_weights.w_U = optimal_weights_array(3);
optimal_weights.w_N = optimal_weights_array(4);
optimal_weights.w_S = optimal_weights_array(5);
optimal_weights.w_T = optimal_weights_array(6);
optimal_weights.w_F = optimal_weights_array(7);

%% Run comparison between original and optimized weights
fprintf('\nRunning comparison between original and optimized weights...\n');

% Run with original weights
original_weights_array = [original_weights.w_P, original_weights.w_K, ...
                         original_weights.w_U, original_weights.w_N, ...
                         original_weights.w_S, original_weights.w_T, original_weights.w_F];

original_fitness = evaluateWeightsComplete(original_weights_array);
                               
fprintf('\nOriginal weights fitness: %.4f\n', original_fitness);
fprintf('Optimised weights fitness: %.4f\n', final_fitness);
if original_fitness ~= 0
    fprintf('Improvement: %.2f%%\n', (original_fitness - final_fitness)/abs(original_fitness) * 100);
else
    fprintf('Improvement: N/A (original fitness is zero)\n');
end

%% Visualise the optimised path
fprintf('\nGenerating and visualizing path with optimised weights...\n');

% Generate optimized path - using global weights
REWARD_WEIGHTS = optimal_weights_array;
num_steps = 50;
optimized_path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery);
REWARD_WEIGHTS = [];  % Clear global weights

% Visualize optimized path
visualizePath(optimized_path, wind_field, aircraft, evaluation_waypoint, 'Optimised Weights Path');

% Generate path with original weights for comparison
REWARD_WEIGHTS = original_weights_array;
original_path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery);
REWARD_WEIGHTS = [];  % Clear global weights

% Visualize original path
visualizePath(original_path, wind_field, aircraft, evaluation_waypoint, 'Original Weights Path');

%% Define the fitness function for the GA (only optimized weights)
function fitness = evaluateWeights(weights_subset)
    % This function evaluates a subset of weights (only those being optimized)

    global original_weights optimize_indices;

    % Reconstruct complete weights array
    weights_array = [original_weights.w_P, original_weights.w_K, original_weights.w_U, ...
                    original_weights.w_N, original_weights.w_S, original_weights.w_T, original_weights.w_F];

    % Update only the weights being optimized
    weights_array(optimize_indices) = weights_subset;

    % Evaluate using the complete weights array
    fitness = evaluateWeightsComplete(weights_array);
end

% %% Complete fitness evaluation function for loitering optimisation
% function fitness = evaluateWeightsComplete(weights)
%     % This function evaluates the fitness of weights for loitering behavior
%     % Focused on maximizing altitude gain while minimizing energy consumption
%     % Lower fitness value is better
% 
%     global REWARD_WEIGHTS wind_field aircraft initial_state; %evaluation_waypoint excluded
% 
%     % Store weights in global variable for calculate_reward.m to use
%     REWARD_WEIGHTS = weights;
% 
%     % Handle potential errors gracefully
%     try
%         % Generate a path with these weights (which will use global REWARD_WEIGHTS)
%         num_steps = 50;  % Number of steps to simulate for evaluation
%         path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps);
% 
%         % Safety check: Verify path data contains no NaN values
%         if any(isnan([path.states.x, path.states.y, path.states.z])) || ...
%            any(isnan([path.states.V, path.states.thrust]))
%             fprintf('NaN values detected in path - assigning high fitness penalty\n');
%             fitness = 1e6;  % Large penalty
%             REWARD_WEIGHTS = [];  % Clear global weights
%             return;
%         end
% 
%         % ======= ENERGY CONSUMPTION CALCULATION =======
%         % Simple battery consumption model
%         base_power = 10;            % Base power draw (W) - systems, avionics, etc.
%         thrust_coefficient = 0.5;   % Power per unit thrust squared (W/N²)
%         dt = 2;                     % Time per step (seconds)
% 
%         % Calculate total energy used over the path
%         total_energy = 0;
%         for i = 1:length(path.states.thrust)
%             % Power = base + coefficient * thrust²
%             power = base_power + thrust_coefficient * path.states.thrust(i)^2;
%             energy_this_step = power * dt;  % Energy in joules
%             total_energy = total_energy + energy_this_step;
%         end
% 
%         % ======= ALTITUDE PERFORMANCE =======
%         % Calculate net altitude change (positive = gain, negative = loss)
%         altitude_change = path.states.z(end) - initial_state.z;
% 
%         % ======= SAFETY CHECK =======
%         % Verify the drone hasn't dropped below a minimum safe altitude
%         min_acceptable_altitude = mission.alt_floor;  % meters above ground
%         if path.states.z(end) < min_acceptable_altitude
%             fprintf('Unsafe altitude detected: final altitude %.1f m\n', path.states.z(end));
%             fitness = 5e5;  % High penalty
%             REWARD_WEIGHTS = [];  % Clear global weights
%             return;
%         end
% 
%         % ======= FITNESS CALCULATION =======
%         % Simple ratio-based fitness for loitering:
%         % - Higher altitude gain is better (negative term)
%         % - Lower energy consumption is better
%         w_altitude = 5.0;    % Weight for altitude change component
%         w_energy = 1.0;      % Weight for energy consumption component
% 
%         % Convert to appropriate units for balanced scaling
%         altitude_gain_normalized = altitude_change / 100;  % Normalised to 100m units
%         energy_normalized = total_energy / 1000;          % Convert to kJ
% 
%         % Calculate fitness (lower is better)
%         fitness = -w_altitude * altitude_gain_normalized + w_energy * energy_normalized;
% 
%         % Output diagnostic information
%         fprintf('Evaluation: Alt change=%.1fm, Energy=%.1fkJ, Fitness=%.4f\n', ...
%                altitude_change, energy_normalized, fitness);
% 
%         % Clear global weights when done
%         REWARD_WEIGHTS = [];
%     catch e
%         % Handle any errors during evaluation
%         fprintf('Error in fitness evaluation: %s\n', e.message);
%         fitness = 1e6;  % Large penalty for errors
% 
%         % Clear global weights even after errors
%         REWARD_WEIGHTS = [];
%     end
% end

% % %% Complete fitness evaluation function for maximum range optimisation
% function fitness = evaluateWeightsComplete(weights)
%     % This function evaluates fitness for maximum range mode
%     % Focuses on maximizing distance traveled per energy unit
%     % Lower fitness value is better
% 
%     global REWARD_WEIGHTS wind_field aircraft initial_state evaluation_waypoint battery;
% 
%     % Store weights in global variable for calculate_reward.m to use
%     REWARD_WEIGHTS = weights;
% 
%     % Handle potential errors gracefully
%     try
%         % Generate a path with these weights
%         num_steps = 100;  % Increased for better range evaluation
%         path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery);
% 
%         % Safety check for NaN values
%         if any(isnan([path.states.x, path.states.y, path.states.z])) || ...
%            any(isnan([path.states.V, path.states.thrust]))
%             fprintf('NaN values detected in path - assigning high fitness penalty\n');
%             fitness = 1e6;  % Large penalty
%             REWARD_WEIGHTS = [];  % Clear global weights
%             return;
%         end
% 
%         % ======= DISTANCE CALCULATION =======
%         % Get actual displacement toward waypoint (useful progress)
%         initial_waypoint_distance = sqrt((initial_state.x - evaluation_waypoint.x)^2 + ...
%                                         (initial_state.y - evaluation_waypoint.y)^2);
%         final_waypoint_distance = sqrt((path.states.x(end) - evaluation_waypoint.x)^2 + ...
%                                       (path.states.y(end) - evaluation_waypoint.y)^2);
%         waypoint_progress = initial_waypoint_distance - final_waypoint_distance;
% 
%         % % ======= ENERGY CALCULATION =======
%         % % Calculate total energy consumed during flight
%         % total_energy_Wh = 0;
%         % for i = 1:length(path.states.thrust)
%         %     % Calculate power for this step
%         %     power = battery.base_power + battery.thrust_coefficient * path.states.thrust(i);
%         %     power = power / battery.discharge_efficiency;  % Account for discharge efficiency
%         % 
%         %     % Convert to energy (Wh)
%         %     dt = 2;  % seconds
%         %     energy_step = power * (dt/3600);  % Convert to Wh
%         %     total_energy_Wh = total_energy_Wh + energy_step;
%         % end
%         % 
%         % % ======= SAFETY CHECKS =======
%         % min_acceptable_altitude = 200;  % meters - using direct value rather than undefined variable
%         % if path.states.z(end) < min_acceptable_altitude
%         %     fprintf('Unsafe altitude detected: final altitude %.1f m\n', path.states.z(end));
%         %     fitness = 5e5;  % High penalty
%         %     REWARD_WEIGHTS = [];  % Clear global weights
%         %     return;
%         % end
% 
%         % ======= FITNESS CALCULATION =======
%         % Calculate altitude gain component (corrected 'start' to '1')
%         alt_gained = path.states.z(end) - path.states.z(1);
% 
%         % Convert altitude to potential distance using 10:1 glide ratio
%         potential_distance = alt_gained * 10;
% 
%         % Total effective distance (actual progress + potential from altitude)
%         effective_distance = waypoint_progress + potential_distance;
% 
%         % % Calculate efficiency metric if energy was consumed
%         % if total_energy_Wh > 0
%         %     distance_per_energy = effective_distance / total_energy_Wh;
%         % else
%         %     distance_per_energy = 0;
%         % end
% 
%         % Final fitness (negative because we want to maximize)
%         % Include both raw distance and efficiency
%         fitness = -(effective_distance); % * 0.7 + distance_per_energy * 100 * 0.3);
% 
%         % % Output diagnostic information
%         % fprintf('Evaluation: Progress=%.2fm, Alt diff=%.2fm, Energy=%.2fWh, Dist/Energy=%.2fm/Wh, Fitness=%.4f\n', ...
%         %        waypoint_progress, alt_gained, total_energy_Wh, distance_per_energy, fitness);
% 
%         fprintf('Evaluation: Progress=%.2fm, Alt diff=%.2fm, Fitness=%.4f\n', ...
%         waypoint_progress, alt_gained, fitness);
% 
%         % Clear global weights when done
%         REWARD_WEIGHTS = [];
%     catch e
%         % Handle any errors during evaluation
%         fprintf('Error in fitness evaluation: %s\n', e.message);
%         fitness = 1e6;  % Large penalty for errors
%         REWARD_WEIGHTS = [];  % Clear global weights
%     end
% end


function fitness = evaluateWeightsComplete(weights)
    % This function evaluates fitness for maximum range and energy optimization
    % Focuses on maximizing distance traveled while maintaining energy reserves
    % Lower fitness value is better

    global REWARD_WEIGHTS wind_field aircraft initial_state evaluation_waypoint battery mission;

    % Store weights in global variable for calculate_reward.m to use
    REWARD_WEIGHTS = weights;

    % Handle potential errors gracefully
    try
        %% Generate path with current weights
        num_steps = 100;  % Increased for better range evaluation
        path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery);

        % Safety check for NaN values
        if any(isnan([path.states.x, path.states.y, path.states.z])) || ...
           any(isnan([path.states.V, path.states.thrust]))
            fprintf('NaN values detected in path - assigning high fitness penalty\n');
            fitness = 1e6;  % Large penalty
            REWARD_WEIGHTS = [];  % Clear global weights
            return;
        end

        %% Waypoint Progress Term
        initial_waypoint_distance = sqrt((initial_state.x - evaluation_waypoint.x)^2 + ...
                                        (initial_state.y - evaluation_waypoint.y)^2);
        final_waypoint_distance = sqrt((path.states.x(end) - evaluation_waypoint.x)^2 + ...
                                      (path.states.y(end) - evaluation_waypoint.y)^2);
        waypoint_progress = initial_waypoint_distance - final_waypoint_distance;
        
        % Normalize waypoint progress (0-1 scale)
        if initial_waypoint_distance > 0
            D_progress = waypoint_progress / initial_waypoint_distance;
        else
            D_progress = 0;
        end

        %% Energy Terms
        % Battery Remaining Term
        if isfield(path.states, 'battery_per')
            fprintf('battery GOOD');
            E_remaining = path.states.battery_per(end) / 100; % Assuming percentage is 0-100 scale
            
        else
            fprintf('baddd no battery ahh!!!!!!!');
        end
        
        % Potential Energy Term (altitude)
        % h_min = mission.alt_floor; % Minimum safe altitude
        h_min = 150;
        fprintf('mission alt floor = %.f ', h_min);
        h_max = 800; % Maximum expected altitude - adjust based on your scenario
        E_potential = max(0, min(1, (path.states.z(end) - h_min) / (h_max - h_min)));
        
        % Kinetic Energy Term (airspeed)
        V_min = aircraft.min_speed;
        V_max = aircraft.max_speed;
        E_kinetic = max(0, min(1, (path.states.V(end) - V_min) / (V_max - V_min)));

        %% Penalty Terms
        % Stall Penalty (minimum speed violation)
        V_min_observed = min(path.states.V_rel);
        if V_min_observed < V_min
            P_stall = 10000 * (V_min - V_min_observed)^2;
        else
            P_stall = 0;
        end
        
        % Altitude Floor Penalty
        if h_min > 0
            fprintf('mission floor GOOD');
            alt_min_observed = min(path.states.z);
            if alt_min_observed < h_min
                P_altitude = 10000 * (h_min - alt_min_observed)^2;
            else
                P_altitude = 0;
            end
        else
            fprintf('mission floor not being read');
        end

        %% Combined Fitness Calculation
        % Weights for each component
        alpha = 30;  % Waypoint progress weight
        beta = 40;   % Battery remaining weight
        gamma = 30;  % Potential energy weight
        delta = 5;  % Kinetic energy weight
        
        % Calculate fitness (negative because we want to maximize the composite term)
        fitness = -(alpha * D_progress + beta * E_remaining + gamma * E_potential + delta * E_kinetic) + P_stall + P_altitude;

        % Output diagnostic information
        fprintf('Eval: Prog=%.2f, Batt=%.1f%%, Alt=%.1fm, Speed=%.1fm/s, Penalties=%.1f, Fitness=%.4f\n', ...
               D_progress, E_remaining * 100, path.states.z(end), path.states.V(end), P_stall + P_altitude, fitness);

        % Clear global weights when done
        REWARD_WEIGHTS = [];
    catch e
        % Handle any errors during evaluation
        fprintf('Error in fitness evaluation: %s\n', e.message);
        fitness = 1e6;  % Large penalty for errors
        REWARD_WEIGHTS = [];  % Clear global weights
    end
end


%% Visualization function
function visualizePath(path, wind_field, aircraft, waypoint, title_text)
    figure;
    
    % Plot wind field
    quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
        wind_field.Wx, wind_field.Wy, wind_field.U, ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.01);  % Light grey color
    hold on
    
    % Plot path with velocity-based coloring
    velocities = path.states.V;
    normalized_velocities = (velocities - aircraft.min_speed) / ...
        (aircraft.max_speed - aircraft.min_speed);
    colors = jet(100);  % Create colormap
    
    % Plot path segments
    for i = 1:length(velocities)-1
        color_idx = max(1, min(100, round(normalized_velocities(i) * 99) + 1));
        plot3([path.states.x(i), path.states.x(i+1)], ...
              [path.states.y(i), path.states.y(i+1)], ...
              [path.states.z(i), path.states.z(i+1)], ...
              'Color', colors(color_idx,:), 'LineWidth', 2);
    end
    
    % Plot waypoint
    plot3(waypoint.x, waypoint.y, waypoint.z, 'r*', ...
        'MarkerSize', 10, ...
        'LineWidth', 2);
    
    % Add colorbar
    c = colorbar;
    ylabel(c, 'Velocity (m/s)')
    % Set colorbar ticks to actual velocity values
    cticks = linspace(0, 1, 5);
    cticklabels = linspace(aircraft.min_speed, aircraft.max_speed, 5);
    c.Ticks = cticks;
    c.TickLabels = arrayfun(@(x) sprintf('%.1f', x), cticklabels, 'UniformOutput', false);
    
    % Mark start point
    plot3(path.states.x(1), path.states.y(1), path.states.z(1), ...
        'go', 'MarkerSize', 10, 'LineWidth', 2)
    
    % Format plot
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    title(title_text)
    grid on
    view(3)
    axis equal
    
    % Calculate energy consumption for statistics
    base_power = 10;     % Power for avionics/systems at idle (W)
    thrust_coefficient = 0.5;  % Power per unit thrust squared (W/N²)
    dt = 2;  % seconds per step
    total_energy = 0;
    for i = 1:length(path.states.thrust)
        power = base_power + thrust_coefficient * path.states.thrust(i)^2;
        total_energy = total_energy + (power * dt);
    end
    
    % Add additional statistics
    thrust_levels = path.states.thrust;
    annotation('textbox', [0.15, 0.05, 0.3, 0.15], ...
        'String', {sprintf('Avg Thrust: %.1f', mean(thrust_levels)), ...
                  sprintf('Zero Thrust: %.1f%%', sum(thrust_levels == 0)/length(thrust_levels)*100), ...
                  sprintf('Final Alt: %.1f m', path.states.z(end)), ...
                  sprintf('Total Energy: %.1f kJ', total_energy/1000), ...
                  sprintf('Dist to WP: %.1f m', sqrt((path.states.x(end) - waypoint.x)^2 + ...
                                                 (path.states.y(end) - waypoint.y)^2 + ...
                                                 (path.states.z(end) - waypoint.z)^2))}, ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'white');
end

%% Inline if function (ternary operator substitute)
function result = iif(condition, true_value, false_value)
    if condition
        result = true_value;
    else
        result = false_value;
    end
end



%% Function to collect GA data
function [state, options, optchanged] = saveGAData(options, state, flag)
    global ga_history original_weights optimize_indices date_folder optimization_id;
    optchanged = false;
    
    % Only process at iteration events
    if strcmp(flag, 'iter')
        % Store generation number
        ga_history.generation(end+1) = state.Generation;
        
        % Store best fitness
        ga_history.best_fitness(end+1) = state.Best(end);
        
        % Store best weights from this generation
        if ~isempty(state.Population)
            [~, best_idx] = min(state.Score);
            best_subset = state.Population(best_idx(1),:);
            
            % Reconstruct full weights array
            best_weights_array = [original_weights.w_P, original_weights.w_K, original_weights.w_U, ...
                                 original_weights.w_N, original_weights.w_S, original_weights.w_T, original_weights.w_F];
            
            % Update optimized weights
            best_weights_array(optimize_indices) = best_subset;
            
            % Store
            ga_history.best_weights(end+1,:) = best_weights_array;
        end
    elseif strcmp(flag, 'done')
        % Save final data
        filename = fullfile(date_folder, [optimization_id '.mat']);
        save(filename, 'ga_history');
        fprintf('Optimization history saved to: %s\n', filename);
    end
end