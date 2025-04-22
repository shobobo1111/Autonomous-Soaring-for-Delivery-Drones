% function [next_state] = flight_dynamics(current_state, control_input, dt, wind, aircraft)
%     % FLIGHT_DYNAMICS Propagates aircraft state using 3-DOF equations from ALFA paper
%     %
%     % Inputs:
%     %   current_state: Structure containing:
%     %       - x, y, z: Position (m)
%     %       - V: Airspeed (m/s)
%     %       - gamma: PITCH angle (rad)
%     %       - phi: ROLL/Bank angle (rad)
%     %       - psi: YAW/Heading angle (rad)
%     %
%     %   control_input: Structure containing:
%     %       - phi_command: Commanded bank angle (rad)
%     %       - thrust: Thrust command (N)
%     %
%     %   dt: Time step (s)
%     %   wind: Structure containing:
%     %       - Wx, Wy: Horizontal wind components (m/s)
%     %       - U: Vertical wind component (m/s)
%     %
%     %   aircraft: Structure containing:
%     %       - mass: Aircraft mass (kg)
%     %       - wing_area: Wing area (m^2)
%     %       - CD0: Zero-lift drag coefficient
%     %       - k: Induced drag factor
%     %
%     % Outputs:
%     %   next_state: Updated aircraft state after dt seconds
% 
%     % Extract current states
%     x = current_state.x;
%     y = current_state.y;
%     z = current_state.z;
%     V = current_state.V;
%     gamma = current_state.gamma;
%     phi = current_state.phi;
%     psi = current_state.psi;
% 
%     % Physical constants
%     g = 9.81;
%     rho = 1.225;
% 
%     % Calculate lift coefficient (assuming steady flight)
%     % CL = 2 * aircraft.mass * g * cos(phi) / (rho * V^2 * aircraft.wing_area);
%     CL = 0.001;
% 
%     % Calculate drag coefficient
%     CD = aircraft.CD0 + aircraft.k * CL^2;
% 
%     % Calculate aerodynamic forces
%     q = 0.5 * rho * V^2 * aircraft.wing_area;
%     % L = q * CL;
%     D = q * CD;
% 
%     % U_v = L - aircraft.mass; 
% 
%     % Flight state parameters
%     nominal_glide_speed = 12;  % m/s
%     powered_speed = 15;        % m/s
%     base_sink_rate = -1;     % m/s
% 
%     % Update velocity based on flight state
%     if control_input.thrust == 0
%         % Unpowered - velocity affected by pitch
%         % Positive gamma (pitch up) reduces speed
%         % Negative gamma (pitch down) increases speed
%         speed_change = -2 * sin(gamma);  % Simple linear relationship
%         V = V + speed_change * dt;
% 
%         % Limit velocity to realistic range for glider
%         V = max(8, min(20, V));  % Between 8-20 m/s
% 
%         % Sink rate might increase at very low or high speeds
%         speed_factor = abs(V - nominal_glide_speed) / nominal_glide_speed;
%         actual_sink_rate = base_sink_rate * (1 + speed_factor);
%     else
%         % Powered flight - maintain constant speed
%         V = powered_speed;
%         actual_sink_rate = 0;
%     end
% 
% 
%     % State derivatives
%     turn_rate = (g * tan(phi)) / V;
%     dpsi = turn_rate;
% 
%     % Simple pitch rate control (proportional control)
%     gamma_error = control_input.gamma_command - gamma;
%     gamma_rate = 1.0 * gamma_error;  % 1.0 is a control gain
% 
%     % Position derivatives using commanded gamma
%     dx = V * cos(gamma) * cos(psi) + wind.Wx;
%     dy = V * cos(gamma) * sin(psi) + wind.Wy;
%     dz = V * sin(gamma) + actual_sink_rate + wind.U;
% 
%     % % Velocity and flight path angle derivatives
%     % dV = (control_input.thrust * cos(0) - D - aircraft.mass * g * sin(gamma)) / aircraft.mass;
% 
%     % Update states
%     next_state.x = x + dx * dt;
%     next_state.y = y + dy * dt;
%     next_state.z = z + dz * dt;
%     next_state.V = V; % + dV * dt;
%     next_state.psi = psi + dpsi * dt;
%     next_state.gamma = gamma + gamma_rate * dt;  % Update gamma based on control
%     next_state.phi = control_input.phi_command;
% end

% function [next_state] = flight_dynamics(current_state, control_input, dt, wind, aircraft)
% % FLIGHT_DYNAMICS Simulates 3-DOF flight dynamics for a fixed-wing UAV
% % Implements both powered and unpowered (gliding) flight modes
% % Based on fundamental aerodynamic principles and energy equations
% 
% % Extract current states
% x = current_state.x;          % East position (m)
% y = current_state.y;          % North position (m)
% z = current_state.z;          % Altitude (m)
% V = current_state.V;          % Airspeed (m/s)
% gamma = current_state.gamma;   % Pitch angle (rad)
% phi = current_state.phi;      % Bank angle (rad)
% psi = current_state.psi;      % Heading angle (rad)
% 
% % Physical constants
% g = 9.81;    % Gravitational acceleration (m/s^2)
% rho = 1.225; % Air density at sea level (kg/m^3)
% 
% % Calculate aerodynamic forces
% % Dynamic pressure q = 1/2 * rho * V^2
% q = 0.5 * rho * V^2 * aircraft.wing_area;
% 
% % Simplified lift coefficient (constant CL assumption)
% % TODO: Implement angle of attack dependence: CL = CL0 + CLα * alpha
% CL = aircraft.CL0;
% 
% % Total drag coefficient (parasitic + induced drag)
% CD = aircraft.CD0 + aircraft.k * CL^2;
% 
% % Calculate drag force
% D = q * CD;
% 
% % Calculate weight
% W = aircraft.mass * g;
% 
% % Calculate minimum sink rate using energy equations
% % V_unpowered_glide: Velocity for minimum sink rate in still air
% % Derived from L = W and minimum power required
% V_unpowered_glide = sqrt((2*W/rho*aircraft.wing_area)*CL);
% V_sink_min = V_unpowered_glide * CD/CL;
% 
% % Current sink rate at actual velocity
% V_sink = V * CD/CL;
% 
% % Debug output
% % fprintf('CL: %f CD: %f Drag: %f Weight: %f\n', CL, CD, D, W);
% % fprintf('V_sink_min: %f\n',V_sink_min);
% % fprintf('V_sink: %f\n',V_sink);
% 
% % % Calculate velocity change (from ALFA paper)
% % dV = (control_input.thrust - D - aircraft.mass * g * sin(gamma)) / aircraft.mass;
% % V = V + dV * dt;
% % 
% % % Ensure velocity stays within limits
% % V = max(aircraft.min_speed, min(aircraft.max_speed, V));
% 
% if control_input.thrust == 0
%     % UNPOWERED FLIGHT (GLIDING)
%     % Simple speed management based on pitch
%     % Negative pitch increases speed, positive pitch decreases
%     % Factor of -2 is a simplified control gain
%     speed_change = -2 * sin(gamma);
%     V = V + speed_change * dt;
%     V = max(aircraft.min_speed, min(aircraft.max_speed, V));
% 
%     % Sink rate varies with deviation from best glide speed
%     speed_factor = abs(V - aircraft.cruise_speed) / aircraft.cruise_speed;
%     actual_sink_rate = V_sink_min * (1 + speed_factor);
% 
% else
%     % POWERED FLIGHT
%     % Calculate excess power available for climb
%     % P = F*V, excess power = thrust power - drag power
%     excess_power = (control_input.thrust * V - D * V) / aircraft.mass;
% 
%     % Convert excess power to vertical velocity
%     % P = m*g*V_vertical -> V_vertical = P/(m*g)
%     actual_sink_rate = excess_power / g; % Negative for climb
% 
%     % Manage airspeed based on thrust/drag balance
%     V = powered_speed_management(V, control_input.thrust, D, aircraft);
% end
% 
% % Calculate turn dynamics
% % From handbook p.3-13: turn rate = g*tan(bank)/velocity
% turn_rate = (g * tan(phi)) / V;
% dpsi = turn_rate;
% 
% % Basic pitch control (proportional control)
% gamma_error = control_input.gamma_command - gamma;
% gamma_rate = 1.0 * gamma_error;  % 1.0 is control gain
% 
% % Update position
% % Include wind effects in ground-relative motion
% dx = V * cos(gamma) * cos(psi) + wind.Wx;
% dy = V * cos(gamma) * sin(psi) + wind.Wy;
% dz = V * sin(gamma) + V_sink + wind.U;
% 
% % Update state vector
% next_state.x = x + dx * dt;
% next_state.y = y + dy * dt;
% next_state.z = max(0, z + dz * dt);  % Prevent negative altitude
% next_state.V = V;
% next_state.psi = psi + dpsi * dt;
% next_state.gamma = gamma + gamma_rate * dt;
% next_state.phi = control_input.phi_command;
% 
% % Store debug info for analysis
% next_state.debug.L = q * CL;
% next_state.debug.D = D;
% next_state.debug.sink_rate = actual_sink_rate;
% end
% 
% function V_new = powered_speed_management(V, thrust, drag, aircraft)
% % Simplified speed management for powered flight
% % Assumes constant acceleration/deceleration rates
% % TODO: Implement proper acceleration based on excess thrust
% if thrust > drag
%     V_new = min(V + 0.1, aircraft.max_speed);  % 0.1 m/s^2 acceleration
% elseif thrust < drag
%     V_new = max(V - 0.1, aircraft.min_speed);  % 0.1 m/s^2 deceleration
% else
%     V_new = V;  % Maintain current speed
% end
% end

% function [next_state] = flight_dynamics(current_state, control_input, dt, wind, aircraft)
% % FLIGHT_DYNAMICS Simulates 3-DOF flight dynamics for a fixed-wing UAV
% % Based on ALFA paper equations and fundamental aerodynamic principles
% 
% % Extract current states
% x = current_state.x;          % East position (m)
% y = current_state.y;          % North position (m)
% z = current_state.z;          % Altitude (m)
% V = current_state.V;          % Airspeed (m/s)
% gamma = current_state.gamma;   % Pitch angle (rad)
% phi = current_state.phi;      % Bank angle (rad)
% psi = current_state.psi;      % Heading angle (rad)
% 
% % Calculate aircraft velocity components
% Vx = V * cos(gamma) * cos(psi);
% Vy = V * cos(gamma) * sin(psi);
% Vz = V * sin(gamma);
% 
% % Get relative wind
% Vx_rel = Vx - wind.Wx;
% Vy_rel = Vy - wind.Wy;
% Vz_rel = Vz - wind.U;
% 
% % Calculate true angle of attack
% V_horiz = max(sqrt(Vx_rel^2 + Vy_rel^2), 0.1);  % Prevent divide by zero
% alpha = gamma - atan2(Vz_rel, V_horiz);
% 
% % Limit alpha to realistic range
% alpha_max = 15 * pi/180;  % ~15 degrees, typical stall angle
% alpha = min(max(alpha, -alpha_max), alpha_max);
% 
% % Using true alpha for aerodynamic calculations
% CL = aircraft.CL0 + aircraft.CL_alpha * alpha;
% CD = aircraft.CD0 + aircraft.k * CL^2;
% 
% 
% % Physical constants
% g = 9.81;    % Gravitational acceleration (m/s^2)
% rho = 1.225; % Air density at sea level (kg/m^3)
% 
% % Calculate aerodynamic forces
% q = 0.5 * rho * V^2 * aircraft.wing_area;
% D = q * CD;
% 
% % Calculate sink rate based on drag polar
% % V_sink = V * CD/CL;
% % V_sink = V_sink * 1/cos(phi); % this increases the sink rate with bank angle
% V_sink = V * sin(gamma) + (g/V) * cos(gamma); % Basic energy balance
% if phi ~= 0
%     V_sink = V_sink * 1/cos(phi); % Bank angle effect
% end
% 
% % fprintf('(%.3f , %.3f)',1/cos(phi),phi);
% % fprintf('V_sink: %f, V: %f', V_sink, V)
% 
% % Calculate velocity change (ALFA paper equation 1)
% dV = (control_input.thrust - D - aircraft.mass * g * sin(gamma)) / aircraft.mass;
% V = V + dV * dt;
% V = max(aircraft.min_speed, min(aircraft.max_speed, V));
% 
% % Calculate turn dynamics
% turn_rate = (g * tan(phi)) / V;
% dpsi = turn_rate;
% 
% % Basic pitch control
% gamma_error = control_input.gamma_command - gamma;
% gamma_rate = 1.0 * gamma_error;
% % lift_force = q * CL;
% % gamma_rate = (1/(aircraft.mass * V)) * (lift_force * cos(phi) - aircraft.mass * g * cos(gamma));
% 
% % Position derivatives including wind and sink rate
% dx = V * cos(gamma) * cos(psi) + wind.Wx;
% dy = V * cos(gamma) * sin(psi) + wind.Wy;
% dz = V * sin(gamma) - V_sink + wind.U;
% 
% % Update state vector
% next_state.x = x + dx * dt;
% next_state.y = y + dy * dt;
% next_state.z = max(0, z + dz * dt);  % Prevent negative altitude
% next_state.V = V;
% next_state.psi = psi + dpsi * dt;
% next_state.gamma = gamma + gamma_rate * dt;
% next_state.phi = control_input.phi_command;
% 
% % Store debug info
% next_state.debug.L = q * CL;
% next_state.debug.D = D;
% next_state.debug.sink_rate = V_sink;
% end

function [next_state] = flight_dynamics(current_state, control_input, dt, wind, aircraft)
% FLIGHT_DYNAMICS Simulates 3-DOF flight dynamics for a fixed-wing UAV
% Follows equations from Powers et al. (2020) ALFA paper

% Extract current states
x = current_state.x;      % East position (m)
y = current_state.y;      % North position (m)
z = current_state.z;      % Altitude (m)
V = current_state.V;      % Airspeed (m/s)
% V_rel = current_state.V_rel;      % Airspeed (m/s)
gamma = current_state.gamma;  % Flight path angle (rad)
phi = current_state.phi;     % Bank angle (rad)
psi = current_state.psi;     % Heading angle (rad)

% Calculate aircraft velocity components
Vx = V * cos(gamma) * cos(psi);
Vy = V * cos(gamma) * sin(psi);
Vz = V * sin(gamma);

% Get relative wind
Vx_rel = Vx - wind.Wx;
Vy_rel = Vy - wind.Wy;
Vz_rel = Vz - wind.U;

% Calculate true angle of attack (α)
V_horiz = max(sqrt(Vx_rel^2 + Vy_rel^2), 0.1);  % Prevent divide by zero

% alpha = gamma - atan2(Vz_rel, V_horiz);
% 
% % Limit alpha to realistic range
% alpha_max = 15 * pi/180;  % ~15 degrees, typical stall angle
% alpha = min(max(alpha, -alpha_max), alpha_max);
% 
% % Calculate aerodynamic coefficients
% CL = aircraft.CL0 + aircraft.CL_alpha * alpha;
% CD = aircraft.CD0 + aircraft.k * CL^2;

% Calculate true angle of attack (α) - Don't cap this one
alpha_true = gamma - atan2(Vz_rel, V_horiz);

% Calculate aerodynamic coefficients with stall modeling
% Reference: Based on typical airfoil behavior described in glider handbook.pdf
% Note: These values should be tuned based on actual aircraft data
alpha_stall = 15 * pi/180;  % Stall angle in radians
alpha_zero_lift = -2 * pi/180;  % Alpha at zero lift

% Pre-stall regime
if abs(alpha_true) <= alpha_stall
    CL = aircraft.CL0 + aircraft.CL_alpha * alpha_true;
    CD = aircraft.CD0 + aircraft.k * CL^2;
else
    % Post-stall regime
    % Smooth transition to post-stall behavior
    % Using a simplified model based on common airfoil characteristics
    alpha_excess = abs(alpha_true) - alpha_stall;
    stall_factor = exp(-alpha_excess/(5*pi/180));  % Exponential decay after stall
    
    % CL drops off after stall but doesn't go to zero immediately
    CL_stall = aircraft.CL0 + aircraft.CL_alpha * alpha_stall;
    CL = sign(alpha_true) * CL_stall * stall_factor;
    
    % CD increases significantly in stall
    CD = aircraft.CD0 + aircraft.k * CL^2 + 2.0 * (1 - stall_factor);
end


% Physical constants
g = 9.81;  % Gravitational acceleration (m/s^2)
rho = 1.225;  % Air density at sea level (kg/m^3)

% Calculate aerodynamic forces
q = 0.5 * rho * V^2 * aircraft.wing_area;
L = q * CL;  % Calculate lift force
% fprintf('Lift before %.3f  ',L);
% L = L*cos(phi); % adjusts Lift when banking to find lift upwards component  
% fprintf('Lift after bank reduction: %.3f \n',L);
D = q * CD;
W = aircraft.mass * g;

% % Apply bank-dependent lift correction
% bank_correction = cos(phi);  % Simple lift reduction based on bank angle
% L = L * bank_correction;  % Reduce lift based on bank angle
% 
% if L < W && control_input.thrust == 0
%     % Force nose down proportional to lift deficit 
%     gamma = -0.1 * (1 - L/W);  % Simple linear relationship
% end

% Calculate velocity change (From Powers et al. Eq 1)
dV = (control_input.thrust*cos(alpha_true) - D - aircraft.mass * g * sin(gamma)) / aircraft.mass;
V = V + dV * dt;
% V = max(aircraft.min_speed, min(aircraft.max_speed, V));

% RE-Calculate aircraft velocity components - giving us V_rel for next
Vx = V * cos(gamma) * cos(psi);
Vy = V * cos(gamma) * sin(psi);
Vz = V * sin(gamma);

% Get relative wind
Vx_rel = Vx - wind.Wx;
Vy_rel = Vy - wind.Wy;
Vz_rel = Vz - wind.U;





% Calculate turn dynamics
turn_rate = (g * tan(phi)) / V;
dpsi = turn_rate;

% Basic pitch control
gamma_error = control_input.gamma_command - gamma;
gamma_rate = 1.0 * gamma_error;
% gamma_rate = (1/(aircraft.mass * V)) * ...
%             ((control_input.thrust * sin(alpha) + ...
%              L )* cos(phi) - ...
%              aircraft.mass * g * cos(gamma));


% Position derivatives - directly from Powers et al. equations
dx = V * cos(gamma) * cos(psi) + wind.Wx;
dy = V * cos(gamma) * sin(psi) + wind.Wy;
dz = V * sin(gamma) + wind.U + (L - W)/(aircraft.mass) * dt;  

% Update state vector
next_state.x = x + dx * dt;
next_state.y = y + dy * dt;
next_state.z = max(0, z + dz * dt);  % Prevent negative altitude
next_state.V = V;
next_state.V_rel = sqrt(Vx_rel^2 + Vy_rel^2 + Vz_rel^2); % used to inform stall in rewards
next_state.psi = psi + dpsi * dt;
next_state.gamma = gamma + gamma_rate * dt;
next_state.phi = control_input.phi_command;


% Store debug info
next_state.debug.L = q * CL;
next_state.debug.D = D;
end


