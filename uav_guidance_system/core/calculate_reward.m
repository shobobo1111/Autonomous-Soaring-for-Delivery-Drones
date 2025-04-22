% function reward = calculate_reward(current_state, next_state, wind, waypoint, elapsed_time)
%     % Constants
%     g = 9.81;
% 
%     % Weights (tune these)
%     w_P = 0.0;    % Weight for potential energy
%     w_U = 0.0;    % Weight for updraft utilisation
%     w_N = 0.5;    % Weight for navigation
% 
%     % Calculate potential energy term (altitude gain)
%     R_P = w_P * g * (next_state.z - current_state.z);
% 
%     % Calculate updraft utiliation term
%     R_U = w_U * 0.5 * abs(wind.U) * wind.U;
% 
%     % Calculate navigation reward with time pressure
%     d_current = sqrt((waypoint.x - current_state.x)^2 + ...
%                     (waypoint.y - current_state.y)^2 + ...
%                     (waypoint.z - current_state.z)^2);
% 
%     d_next = sqrt((waypoint.x - next_state.x)^2 + ...
%                   (waypoint.y - next_state.y)^2 + ...
%                   (waypoint.z - next_state.z)^2);
% 
%     % Time pressure increases navigation reward over time
%     time_factor = waypoint.time_pressure_base + ...
%                  waypoint.time_pressure_rate * elapsed_time;
% 
%     R_N = w_N * (d_current - d_next) * time_factor;
% 
%     % Total reward
%     reward = R_P + R_U + R_N;
% end

% function reward = calculate_reward(current_state, next_state, wind, waypoint, elapsed_time, control_input)
%     global REWARD_WEIGHTS; % Add this line
% 
%     % % Load aircraft configuration
%     aircraft = init_aircraft();
%     mission = init_mission();
% 
%     g = 9.81;
% 
%     % % Set default weights 
%     % w_P = 0.5;    % Weight for potential energy -- 1.2 orginally
%     % w_K = 1.7;    % Weight for kinetic energy -1.5 orginally
%     % w_U = 3.0;    % Weight for updraft utilisation (set to 0 for straight to point flight)
%     % w_N = 0.01;   % Weight for navigation (set to 1 for straight to point flight)
%     % w_S = 5.0;    % Weight for stall prevention
%     % w_T = 0.1;    % Weight for thrust/energy penalty (set to 0 for straight to point flight)
%     % w_F = 1;      % Weight for breaching alt floor
% 
%     % test weights - experimental
%     w_P = 0.2274;    % Weight for potential energy -- 1.2 orginally
%     w_K = 0.8933;    % Weight for kinetic energy -1.5 orginally
%     w_U = 5;    % Weight for updraft utilisation (set to 0 for straight to point flight)
%     w_N = 0.5;   % was good generally at 0.1
%     w_S = 5.0;    % Weight for stall prevention
%     w_T = 0.1;    % Weight for thrust/energy penalty (set to 0 for straight to point flight)
%     w_F = 1;      % Weight for breaching alt floor
% 
% 
%     %new long distance optimised parameters
% 
%     % w_P = 0.2274;
%     % w_K = 0.8933;
%     % w_U = 4.9160;
%     % w_N = 1.0770;
%     % w_S = 5;
%     % w_T = 0.0169;
%     % w_F = 1;
% 
% % % long distance optimised parameters
% %     w_P = 0.5000; % this was 1.2, change back in a bit
% %     w_K = 1.7000; % this was 1.5
% %     w_U = 4.8482;
% %     w_N = 0.0375;
% %     w_S = 5;
% %     w_T = 0.0833;
% %     w_F = 1;
% 
%     % good for long distance (manually tested):
%     % 
%     % w_P = 1.2;    % Weight for potential energy
%     % w_K = 1.5;    % Weight for kinetic energy
%     % w_U = 3.0;    % Weight for updraft utilisation (set to 0 for straight to point flight)
%     % w_N = 0.01;    % Weight for navigation (set to 0 for straight to point flight)
%     % w_S = 5.0;    % Weight for stall prevention
%     % w_T = 0.1;    % Weight for thrust/energy penalty (set to 0 for straight to point flight)
%     % w_F = 1;      % Weight for breaching alt floor
% 
%     % % good for loitering manually tested
%     % w_P = 0.2274;    % Weight for potential energy -- 1.2 orginally
%     % w_K = 0.8933;    % Weight for kinetic energy -1.5 orginally
%     % w_U = 5;    % Weight for updraft utilisation (set to 0 for straight to point flight)
%     % w_N = 0.1;   % granted that waypoint.loiter_within_r = 300
%     % w_S = 5.0;    % Weight for stall prevention
%     % w_T = 0.1;    % Weight for thrust/energy penalty (set to 0 for straight to point flight)
%     % w_F = 1;      % Weight for breaching alt floor
% 
% 
% %     Loitering Optimisation Results: REDO
% %     Optimal weights found:
% %       w_P (potential energy): 1.2000 (fixed)
% %       w_K (kinetic energy): 1.5000 (fixed)
% %       w_U (updraft): 4.7700 
% %       w_N (navigation): 0.1000 - could remove entirely 
% %       w_S (stall prevention): 5.0000 (fixed)
% %       w_T (energy/thrust): 0.1000 
%         % w_P = 1.2;
%         % w_K = 1.5;
%         % w_U = 4.77;
%         % w_N = 0.01;
%         % w_S = 5;
%         % w_T = 0.1;
%         % w_F = 1;
% %     Final fitness value: -2.1694
% 
% 
% 
%     % Override with global weights if provided (for optimization)
%     if ~isempty(REWARD_WEIGHTS)
%         w_P = REWARD_WEIGHTS(1);
%         w_K = REWARD_WEIGHTS(2);
%         w_U = REWARD_WEIGHTS(3);
%         w_N = REWARD_WEIGHTS(4);
%         w_S = REWARD_WEIGHTS(5);
%         w_T = REWARD_WEIGHTS(6);
%         w_F = REWARD_WEIGHTS(7);
%     end
% 
%     alt_floor_w_tol = mission.alt_floor + mission.alt_floor_tol;
%     % 
%     % % creating an exponential penality for breaching the altitide floor
%     % if current_state.z > alt_floor_w_tol
%     %     Dw_P = 1; % dynamic weight inactive till (floor + tol) breached
%     % else 
%     %     Dw_P = (alt_floor_w_tol - current_state.z); 
%     % end
% 
%     % Add altitude floor penalty
%     if next_state.z < alt_floor_w_tol
%         alt_floor_penalty = -10 * (alt_floor_w_tol - next_state.z)^2;
%     % elseif next_state.z < alt_floor_w_tol
%     %     alt_floor_penalty = -5 * (alt_floor_w_tol - next_state.z);
%     else
%         alt_floor_penalty = 0;
%     end
% 
%     R_F = w_F * alt_floor_penalty;
% 
% 
%     % Calculate potential energy term (altitude change)
%     R_P = w_P * g * (next_state.z - current_state.z);
% 
%     % Calculate kinetic energy term (velocity change)
%     R_K = w_K * 0.5 * (next_state.V^2 - current_state.V^2);
% 
%     % Linear stall prevention with soft boundaries
%     stall_margin = 4.0; % m/s buffer above min speed
%     stall_threshold = aircraft.min_speed + stall_margin;
% 
%     if next_state.V_rel < stall_threshold
%         % Linear penalty that increases as speed decreases
%         speed_deficit = stall_threshold - next_state.V_rel;
%         R_S = -w_S * speed_deficit;  % Linear instead of quadratic
%     else
%         R_S = 0;
%     end
% 
%     % Calculate updraft utilisation term
%     R_U = w_U * 0.5 * abs(wind.U) * wind.U;
% 
%     % Calculate navigation reward with time pressure
%     d_current = sqrt((waypoint.x - current_state.x)^2 + ...
%                     (waypoint.y - current_state.y)^2 + ...
%                     (waypoint.z - current_state.z)^2);
%     d_next = sqrt((waypoint.x - next_state.x)^2 + ...
%                   (waypoint.y - next_state.y)^2 + ...
%                   (waypoint.z - next_state.z)^2);
% 
%     given_time = mission.time_to_deliver; % the time allocated for the delivery
%     time_left = given_time - elapsed_time; % given 100s - 20s past already -> you have 80s left
% 
%     % the predicted time it will take to get to the waypoint at cruise
%     % added buffer to warn of lateness before it happens
%     predicted_time = d_current / aircraft.cruise_speed; % 150m left to go / 15m/s -> 10s predicted time to get there
%     p_time_w_buffer = predicted_time + (predicted_time * 0.5); % now 15s to get there cause life
% 
%     % lateness is how late you would arrive at cruise (plus buffer)
%     lateness = p_time_w_buffer - time_left; % if predicted time is 35 and time left is 30, lateness is 5
% 
%     if lateness > 0
%         time_factor = lateness;
%     else
%         time_factor = 1;
%     end
% 
%     d_current_lin = sqrt((waypoint.x - current_state.x)^2 + ...
%                 (waypoint.y - current_state.y)^2);
% 
%     if d_current_lin > waypoint.loiter_within_r
%         R_N = w_N * (d_current - d_next) * time_factor;
%     else
%         R_N = 0;
%     end
% 
%     % R_N = w_N * (d_current - d_next) * time_factor; % without loitering
% 
%     % dist_closer = (d_current - d_next);
%     % fprintf('lateness %.2f  \ndist closer to waypoint %.2f  \nR_N: %.2f\n\n', lateness, dist_closer, R_N);
% 
%     % Calculate battery consumption penalty based on thrust
%     % if exist('control_input', 'var') && isfield(control_input, 'thrust')
% 
%     % Calculate power consumption in watts based on thrust level
%     % This is a simplified model - you can replace with more accurate data if available
%     % Assuming a quadratic relationship between thrust and power consumption
%     % P = k * T^2, where k is a constant determined by the specific UAV
%     % this is just an estimation to motivate the guidance to reduce power use
% 
%     % Basic battery consumption model parameters
%     base_power = 10;     % Power for avionics/systems at idle (W)
%     thrust_coefficient = 0.5;  % Power per unit thrust squared (W/N²)
% 
%     % Calculate power for this thrust level
%     % power_consumption = base_power + thrust_coefficient * control_input.thrust^2; 
%     % UPDATE, control input thrust will keep its sign, this term will be a
%     % reward not a penality, negative thrust will represent regen
%     power_consumption = base_power + thrust_coefficient * abs(control_input.thrust) * control_input.thrust;
% 
%     % Calculate energy used in this time step (2 seconds as per your simulation)
%     dt = 2;  % seconds per step
%     energy_consumed = power_consumption * dt;  % Watt-seconds
% 
%     % Apply penalty based on energy consumption
%     R_T = -w_T * energy_consumed;
% 
%     % else
%     %     R_T = 0; % No penalty if thrust information is unavailable
%     %     fprintf('ahhhhhhhhhh no thrust info');
%     % end
% 
% 
% %     fprintf('Debug Values:\n');
% % fprintf('R_P: %.2f\n', R_P);
% % fprintf('R_K: %.2f\n', R_K);
% % fprintf('R_U: %.2f\n', R_U);
% % fprintf('R_N: %.2f\n', R_N);
% % fprintf('R_S: %.2f\n', R_S);
% % fprintf('R_T: %.2f\n', R_T);
% % fprintf('R_F: %.2f\n', R_F);
% % fprintf('Current Speed: %.2f\n', current_state.V);
% % fprintf('Next Speed: %.2f\n', next_state.V);
% % fprintf('Current Alt: %.2f\n', current_state.z);
% % fprintf('Next Alt: %.2f\n', next_state.z);
% % fprintf('Total Reward: %.2f\n', R_P + R_K + R_U + R_N + R_S + R_T + R_F);
% % fprintf('------------------------\n');
% 
%     % Total reward including stall prevention
%     reward = R_P + R_K + R_U + R_N + R_S + R_T + R_F;
% end
% 


function reward = calculate_reward(current_state, next_state, wind, waypoint, elapsed_time, control_input, operating_mode)
    global REWARD_WEIGHTS; % For optimization

    % Load aircraft configuration
    aircraft = init_aircraft();
    mission = init_mission();

    g = 9.81;

    % Manual wind field boundaries
    x_limits = [-100, 2000];
    y_limits = [-500, 500];
    z_limits = [0, 800];
    

    
    
    % Default to balanced mode if not specified
    if nargin < 7
        operating_mode = 1; % loitering
    end
    
    % Define weight sets for different operating modes
    % Format: [w_P, w_K, w_U, w_N, w_S, w_T, w_F]
    
    % MODE 1: LOITERING
    % Manually tested weights
    loiter_weights = [
        0.2274,  % w_P: Weight for potential energy --- was 0.2274
        0.8933,  % w_K: Weight for kinetic energy --- 0.8933
        5.0,     % w_U: Weight for updraft utilisation
        0.1,     % w_N: Weight for navigation (with loiter radius of 300m)
        5.0,     % w_S: Weight for stall prevention
        0.1,     % w_T: Weight for thrust/energy penalty
        1.0      % w_F: Weight for breaching alt floor
    ];
    
    % Optimisation results for loitering - commented out as alternative
    % loiter_weights = [
    %     1.2,     % w_P: Weight for potential energy (fixed)
    %     1.5,     % w_K: Weight for kinetic energy (fixed)
    %     4.77,    % w_U: Weight for updraft utilisation
    %     0.01,    % w_N: Weight for navigation
    %     5.0,     % w_S: Weight for stall prevention (fixed)
    %     0.1,     % w_T: Weight for thrust/energy penalty
    %     1.0      % w_F: Weight for breaching alt floor
    % ];
    
    % MODE 2: URGENT/STRAIGHT TO POINT
    % No specific weights yet - using default with higher navigation priority
    urgent_weights = [
        0.5,     % w_P: Weight for potential energy
        1.7,     % w_K: Weight for kinetic energy
        0.0,     % w_U: Weight for updraft utilisation (0 for direct flight)
        2.0,     % w_N: Weight for navigation (high for direct flight)
        50.0,     % w_S: Weight for stall prevention
        0.0,     % w_T: Weight for thrust/energy penalty (0 for direct flight)
        1.0      % w_F: Weight for breaching alt floor
    ]; % NOTE: Update later with optimized weights for urgent mode
    
    % MODE 3: LONG DISTANCE

    % % optimised: normalised and weighted terms -> higher p, lower k and battery cap
    % % increased with base power decrease, Final fitness value: -92.4361
    % long_dist_weights = [
    %     1.1450,     % w_P: Weight for potential energy
    %     0.3088,     % w_K: Weight for kinetic energy
    %     4.4524,     % w_U: Weight for updraft utilisation ORGINALLY -> 1.4784
    %     0.2683,    % w_N: Weight for navigation ORGINALLY -> 0.8270
    %     5.0,     % w_S: Weight for stall prevention
    %     0.0169,     % w_T: Weight for thrust/energy penalty
    %     1.0      % w_F: Weight for breaching alt floor
    % ];

    % % optimised: normalised and weighted terms for fitness function
    % long_dist_weights = [
    %     0.1,     % w_P: Weight for potential energy
    %     0.4193,     % w_K: Weight for kinetic energy
    %     4.4784,     % w_U: Weight for updraft utilisation ORGINALLY -> 1.4784
    %     0.28270,    % w_N: Weight for navigation ORGINALLY -> 0.8270
    %     5.0,     % w_S: Weight for stall prevention
    %     0.0169,     % w_T: Weight for thrust/energy penalty
    %     1.0      % w_F: Weight for breaching alt floor
    % ];



    % Manually tested weights
    % long_dist_weights = [
    %     0.5,     % w_P: Weight for potential energy
    %     1.7,     % w_K: Weight for kinetic energy
    %     3.0,     % w_U: Weight for updraft utilisation
    %     0.0375,    % w_N: Weight for navigation
    %     5.0,     % w_S: Weight for stall prevention
    %     0.833,     % w_T: Weight for thrust/energy penalty
    %     1.0      % w_F: Weight for breaching alt floor
    % ];

    % Optimal weights found:
    %   w_P (potential energy): 0.1000 
    %   w_K (kinetic energy): 0.4193 
    %   w_U (updraft): 1.4784 
    %   w_N (navigation): 0.8270 
    %   w_S (stall prevention): 5.0000 (fixed)
    %   w_T (energy/thrust): 0.0169 
    %   w_F (alt floor breach): 1.0000 (fixed)
    
    % Optimisation results for long distance - commented out as alternative
    % long_dist_weights = [
    %     0.5000,  % w_P: Weight for potential energy
    %     1.7000,  % w_K: Weight for kinetic energy
    %     4.8482,  % w_U: Weight for updraft utilisation
    %     0.0375,  % w_N: Weight for navigation
    %     5.0,     % w_S: Weight for stall prevention
    %     0.0833,  % w_T: Weight for thrust/energy penalty
    %     1.0      % w_F: Weight for breaching alt floor
    % ];
    
    % MODE 4: BALANCED ----------- NOT USING
    % These are the "test weights - experimental" from your code
    balanced_weights = [
        0.2274,  % w_P: Weight for potential energy 
        0.8933,  % w_K: Weight for kinetic energy 
        5.0,     % w_U: Weight for updraft utilisation
        0.1,     % w_N: Weight for navigation 
        10.0,     % w_S: Weight for stall prevention
        0.1,     % w_T: Weight for thrust/energy penalty
        1.0      % w_F: Weight for breaching alt floor
    ];
    
    % MODE 5: ENERGY SAVING 
    energy_saving_weights = [
        0.2274,  % w_P: Weight for potential energy
        0.8933,  % w_K: Weight for kinetic energy
        5.0,     % w_U: Weight for updraft utilisation (high to maximize energy saving)
        0.1,     % w_N: Weight for navigation (low for energy saving)
        5.0,     % w_S: Weight for stall prevention
        0.1,     % w_T: Weight for thrust/energy penalty (higher penalty for energy use)
        1.0      % w_F: Weight for breaching alt floor
    ]; % NOTE: Update later with optimized weights for energy saving mode
    
    % MODE 6: REGENERATION
    regen_weights = [
        0.82274,  % w_P: Weight for potential energy
        0.8933,  % w_K: Weight for kinetic energy
        5.0,     % w_U: Weight for updraft utilisation
        0.1,     % w_N: Weight for navigation
        5.0,     % w_S: Weight for stall prevention
        0.1,    % w_T: Weight for thrust/energy penalty 
        1.0      % w_F: Weight for breaching alt floor
    ]; % NOTE: Update later with optimized weights for regeneration mode
    
    % Select weights based on operating mode
    switch operating_mode
        case 1  % Loiter
            weights = loiter_weights;
        case 2  % Urgent
            weights = urgent_weights;
        case 3  % Long Distance
            weights = long_dist_weights;
        case 4  % Balanced
            weights = balanced_weights;
        case 5  % Energy Saving
            weights = energy_saving_weights; % gonna rename to 
        case 6  % Regeneration
            weights = regen_weights;
        otherwise
            weights = energy_saving_weights;  % Default to balanced
    end
    
    % Extract individual weights
    w_P = weights(1);
    w_K = weights(2);
    w_U = weights(3);
    w_N = weights(4);
    w_S = weights(5);
    w_T = weights(6);
    w_F = weights(7);
    
    % Override with global weights if provided (for optimization)
    if ~isempty(REWARD_WEIGHTS)
        w_P = REWARD_WEIGHTS(1);
        w_K = REWARD_WEIGHTS(2);
        w_U = REWARD_WEIGHTS(3);
        w_N = REWARD_WEIGHTS(4);
        w_S = REWARD_WEIGHTS(5);
        w_T = REWARD_WEIGHTS(6);
        w_F = REWARD_WEIGHTS(7);
    end

    alt_floor_w_tol = mission.alt_floor + mission.alt_floor_tol;

    % Add altitude floor penalty
    if next_state.z < alt_floor_w_tol
        alt_floor_penalty = -10 * (alt_floor_w_tol - next_state.z)^2;
    else
        alt_floor_penalty = 0;
    end

    R_F = w_F * alt_floor_penalty;

    % Calculate potential energy term (altitude change)
    R_P = w_P * g * (next_state.z - current_state.z);
    
    % Calculate kinetic energy term (velocity change)
    R_K = w_K * 0.5 * (next_state.V^2 - current_state.V^2);
    
    % Linear stall prevention with soft boundaries
    stall_margin = 4.0; % m/s buffer above min speed
    stall_threshold = aircraft.min_speed + stall_margin;
    
    if next_state.V_rel < stall_threshold
        % Linear penalty that increases as speed decreases
        speed_deficit = stall_threshold - next_state.V_rel;
        R_S = -w_S * speed_deficit;  % Linear instead of quadratic
    else
        R_S = 0;
    end
    
    % Calculate updraft utilisation term
    % updraft_v = sqrt(abs(wind.U) * wind.U);
    
    R_U = w_U * 0.5 * wind.U; % this used be not sqrt, causing exponential rewards to go into hella strong updrafts
    % updraft_v = sqrt(abs(wind.U) * wind.U);
    % fprintf(' updraft velocity here: %2.f', updraft_v);
    
    % Calculate navigation reward with time pressure
    d_current = sqrt((waypoint.x - current_state.x)^2 + ...
                    (waypoint.y - current_state.y)^2 + ...
                    (waypoint.z - current_state.z)^2);
    d_next = sqrt((waypoint.x - next_state.x)^2 + ...
                  (waypoint.y - next_state.y)^2 + ...
                  (waypoint.z - next_state.z)^2);
                  
    given_time = mission.time_to_deliver; % the time allocated for the delivery
    time_left = given_time - elapsed_time; % given 100s - 20s past already -> you have 80s left
    
    % the predicted time it will take to get to the waypoint at cruise
    % added buffer to warn of lateness before it happens
    predicted_time = d_current / aircraft.cruise_speed; % 150m left to go / 15m/s -> 10s predicted time to get there
    p_time_w_buffer = predicted_time + (predicted_time * 0.5); % now 15s to get there cause life

    % lateness is how late you would arrive at cruise (plus buffer)
    lateness = p_time_w_buffer - time_left; % if predicted time is 35 and time left is 30, lateness is 5

    if lateness > 0
        time_factor = lateness;
    else
        time_factor = 1;
    end
    
    d_current_lin = sqrt((waypoint.x - current_state.x)^2 + ...
                (waypoint.y - current_state.y)^2);

    if d_current_lin > waypoint.loiter_within_r
        R_N = w_N * (d_current - d_next) * time_factor;
    else
        R_N = 0;
    end

    % Calculate battery consumption penalty based on thrust
    % Basic battery consumption model parameters
    base_power = 10;     % Power for avionics/systems at idle (W)
    thrust_coefficient = 0.5;  % Power per unit thrust squared (W/N²)
    
    % Calculate power for this thrust level
    power_consumption = base_power + thrust_coefficient * abs(control_input.thrust) * control_input.thrust;
    
    % Calculate energy used in this time step (2 seconds as per your simulation)
    dt = 2;  % seconds per step
    energy_consumed = power_consumption * dt;  % Watt-seconds
    
    % Apply penalty based on energy consumption
    R_T = -w_T * energy_consumed;


    % Calculate boundary penalty
    boundary_penalty = calculate_boundary_penalty(next_state, x_limits, y_limits, z_limits);
    
    % Add a weight for the boundary penalty - make this high enough to strongly discourage boundary crossing
    w_B = 10.0;  % This value can be adjusted based on testing
    R_B = w_B * boundary_penalty;

    % Total reward including all components
    reward = R_P + R_K + R_U + R_N + R_S + R_T + R_F + R_B;

    % fprintf('Debug Values:\n');
    % fprintf('R_P: %.2f\n', R_P);
    % fprintf('R_K: %.2f\n', R_K);
    % fprintf('R_U: %.2f\n', R_U);
    % fprintf('R_N: %.2f\n', R_N);
    % fprintf('R_S: %.2f\n', R_S);
    % fprintf('R_T: %.2f\n', R_T);
    % fprintf('R_F: %.2f\n', R_F);
    % fprintf('Current Speed: %.2f\n', current_state.V);
    % fprintf('Next Speed: %.2f\n', next_state.V);
    % fprintf('Current Alt: %.2f\n', current_state.z);
    % fprintf('Next Alt: %.2f\n', next_state.z);
    % fprintf('Total Reward: %.2f\n', R_P + R_K + R_U + R_N + R_S + R_T + R_F);
    % fprintf('------------------------\n');
end



function penalty = calculate_boundary_penalty(state, x_limits, y_limits, z_limits)
    % Define safety margin (distance from boundary where penalty starts)
    margin = 50;  % meters
    
    % Initialize penalty
    penalty = 0;
    
    % Check x boundary
    if state.x < x_limits(1) + margin
        % Exponential penalty that grows as aircraft approaches boundary
        penalty = penalty - exp(10 * (1 - (state.x - x_limits(1)) / margin));
    elseif state.x > x_limits(2) - margin
        penalty = penalty - exp(10 * (1 - (x_limits(2) - state.x) / margin));
    end
    
    % Check y boundary
    if state.y < y_limits(1) + margin
        penalty = penalty - exp(10 * (1 - (state.y - y_limits(1)) / margin));
    elseif state.y > y_limits(2) - margin
        penalty = penalty - exp(10 * (1 - (y_limits(2) - state.y) / margin));
    end
    
    % Check z boundary (upper only, lower is already handled by altitude floor)
    if state.z > z_limits(2) - margin
        penalty = penalty - exp(10 * (1 - (z_limits(2) - state.z) / margin));
    end
    
    % return penalty;
end