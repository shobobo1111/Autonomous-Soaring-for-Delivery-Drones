% create_thermals.m
function wind_field = create_thermals(wind_field, thermal_type)
    % Creates different types of thermal patterns in the wind field
    % 
    % Inputs:
    % wind_field - initialized wind field structure
    % thermal_type - 'simple', 'complex', or 'drifting', 'twisting',
    % 'thermal_with_sink', 'sink_drift_complex'
    
    % Default parameters
    params.center = [40, 50];
    params.strength = 6;        % m/s vertical
    params.radius = 30;         % meters
    params.constant_wind = [0, 0, 0];  % [Wx, Wy, Wz]
    
    switch thermal_type
        case 'simple'
            wind_field = create_simple_thermal(wind_field, params);
        case 'complex'
            params.turbulence = 2;    % % variation
            params.noise_scale = 1;    % meter spatial scale
            wind_field = create_complex_thermal(wind_field, params);
        case 'drifting'
            params.drift_direction = [1, 0];  % Drift east
            params.drift_rate = 0.4;    % 0.1m drift per meter height
            wind_field = create_drifting_thermal(wind_field, params);
        case 'twisting'
            params.turbulence = 0.3;       % 30% variation in all directions
            params.noise_scale = 50;       % 50m spatial scale
            params.horiz_strength = 0.5;   % Horizontal wind strength relative to vertical
            wind_field = create_twisting_thermal(wind_field, params);
        case 'thermal_with_sink'
            params.core_radius = 23;       % Core thermal radius (m)
            params.sink_radius = 30;       % Outer sink radius (m)
            params.sink_strength = -1;     % Sink strength relative to thermal (-1 = 20% of thermal strength)
            wind_field = create_thermal_with_sink(wind_field, params);
        case 'sink_drift_complex'
            params.core_radius = 20;
            params.sink_radius = 30;
            params.sink_strength = -0.5;
            params.turbulence = 0.0;
            params.noise_scale = 10;
            params.drift_direction = [1, 0];
            params.drift_rate = 0.4;
            params.horiz_strength = 0;
            wind_field = create_sink_drift_complex(wind_field, params);
        case 'simple_sys'
            % Parameters based on Allen (2005) thermal distribution equations
            params.num_thermals = 25;  % Default, should be calculated from eq 20 when z_i available
            params.rng_seed = 42;     % Fixed random seed for reproducible thermal positions
            params.strength_range = [1, 5];  % Range of thermal strengths in m/s
            params.radius_range = [20, 80];    % Range of thermal radii in m
            % Default winds
            params.constant_wind;  % [Wx, Wy, Wz]
            wind_field = create_simple_thermal_system(wind_field, params);
        case 'drifting_sys'
            % Parameters based on Allen (2005) thermal distribution equations
            params.num_thermals = 40;  % 18 for a 2000x600 area, with a height of 700m and a average raduis of 40
            params.rng_seed = 41;     % Fixed random seed for reproducible thermal positions
            params.strength_range = [1, 5];  % Range of thermal strengths in m/s
            params.radius_range = [20, 60];  % Range of thermal radii in m
            % Drifting parameters
            params.drift_rate = 0.4;  % 0.4m drift per meter height
            % Max height parameter
            params.max_thermal_height = 700;  % Thermals dissipate above this height (m)
            % Use horizontal wind components for drift direction
            params.drift_direction = [params.constant_wind(1), params.constant_wind(2)];
            % Normalize drift direction if not zero
            if norm(params.drift_direction) > 0
                params.drift_direction = params.drift_direction / norm(params.drift_direction);
            end
            wind_field = create_drifting_thermal_system(wind_field, params);
        otherwise
            error('Unknown thermal type');
    end
end

% Local function definitions
function wind_field = create_simple_thermal(wind_field, params)
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                r = sqrt((wind_field.X(j,i,k) - params.center(1))^2 + ...
                        (wind_field.Y(j,i,k) - params.center(2))^2);
                
                wind_field.Wx(j,i,k) = params.constant_wind(1);
                wind_field.Wy(j,i,k) = params.constant_wind(2);
                wind_field.U(j,i,k) = params.strength * ...
                    exp(-(r/params.radius)^2) + params.constant_wind(3);
            end
        end
    end
end

function wind_field = create_complex_thermal(wind_field, params)
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                r = sqrt((wind_field.X(j,i,k) - params.center(1))^2 + ...
                        (wind_field.Y(j,i,k) - params.center(2))^2);
                
                % Add spatial variation
                theta = atan2(wind_field.Y(j,i,k) - params.center(2), ...
                            wind_field.X(j,i,k) - params.center(1));
                radius_variation = params.radius * ...
                    (1 + params.turbulence * sin(theta * 3) * ...
                    sin(wind_field.Z(j,i,k)/params.noise_scale));
                
                % Add strength variations
                strength_variation = params.strength * ...
                    (1 + params.turbulence * ...
                    (sin(wind_field.X(j,i,k)/params.noise_scale) + ...
                     sin(wind_field.Y(j,i,k)/params.noise_scale) + ...
                     sin(wind_field.Z(j,i,k)/params.noise_scale))/3);
                
                wind_field.Wx(j,i,k) = params.constant_wind(1);
                wind_field.Wy(j,i,k) = params.constant_wind(2);
                wind_field.U(j,i,k) = strength_variation * ...
                    exp(-(r/radius_variation)^2) + params.constant_wind(3);
            end
        end
    end
end

function wind_field = create_drifting_thermal(wind_field, params)
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                height = wind_field.Z(j,i,k);
                drift_x = params.drift_direction(1) * height * params.drift_rate;
                drift_y = params.drift_direction(2) * height * params.drift_rate;
                
                r = sqrt((wind_field.X(j,i,k) - (params.center(1) + drift_x))^2 + ...
                        (wind_field.Y(j,i,k) - (params.center(2) + drift_y))^2);
                
                wind_field.Wx(j,i,k) = params.constant_wind(1);
                wind_field.Wy(j,i,k) = params.constant_wind(2);
                wind_field.U(j,i,k) = params.strength * ...
                    exp(-(r/params.radius)^2) + params.constant_wind(3);
            end
        end
    end
end

function wind_field = create_twisting_thermal(wind_field, params)
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                % Base distance from thermal center
                dx = wind_field.X(j,i,k) - params.center(1);
                dy = wind_field.Y(j,i,k) - params.center(2);
                r = sqrt(dx^2 + dy^2);
                
                % Create horizontal variations based on position
                angle = atan2(dy, dx);
                height_factor = wind_field.Z(j,i,k)/params.noise_scale;
                
                % Add rotational component that varies with height
                Wx_rot = -dy/r * sin(height_factor);
                Wy_rot = dx/r * sin(height_factor);
                
                % Add turbulent variations
                turb_x = params.turbulence * sin(dx/params.noise_scale + height_factor);
                turb_y = params.turbulence * cos(dy/params.noise_scale + height_factor);
                
                % Combine components
                wind_field.Wx(j,i,k) = params.horiz_strength * params.strength * ...
                    (Wx_rot + turb_x) * exp(-(r/params.radius)^2) + params.constant_wind(1);
                wind_field.Wy(j,i,k) = params.horiz_strength * params.strength * ...
                    (Wy_rot + turb_y) * exp(-(r/params.radius)^2) + params.constant_wind(2);
                wind_field.U(j,i,k) = params.strength * ...
                    exp(-(r/params.radius)^2) * (1 + 0.2*sin(height_factor)) + params.constant_wind(3);
            end
        end
    end
end

function wind_field = create_thermal_with_sink(wind_field, params)
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                % Distance from thermal center
                r = sqrt((wind_field.X(j,i,k) - params.center(1))^2 + ...
                        (wind_field.Y(j,i,k) - params.center(2))^2);
                
                % Core thermal updraft (Gaussian profile)
                core_updraft = exp(-(r/params.core_radius)^2);
                
                % Surrounding sink (ring of downdraft)
                % Uses difference of two Gaussians to create ring
                sink_profile = exp(-(r/params.sink_radius)^2) - exp(-(r/params.core_radius)^2);
                
                % Combine thermal and sink
                wind_field.Wx(j,i,k) = params.constant_wind(1);
                wind_field.Wy(j,i,k) = params.constant_wind(2);
                wind_field.U(j,i,k) = params.strength * core_updraft + ...
                    params.strength * params.sink_strength * sink_profile + ...
                    params.constant_wind(3);
            end
        end
    end
end

function wind_field = create_sink_drift_complex(wind_field, params)
    % Combines features of complex, drifting, and thermal_with_sink
    
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                % Position and height
                height = wind_field.Z(j,i,k);
                dx = wind_field.X(j,i,k) - params.center(1);
                dy = wind_field.Y(j,i,k) - params.center(2);
                r = sqrt(dx^2 + dy^2);
                
                % Drift position based on altitude
                drift_x = params.drift_direction(1) * height * params.drift_rate;
                drift_y = params.drift_direction(2) * height * params.drift_rate;
                drifted_r = sqrt((dx - drift_x)^2 + (dy - drift_y)^2);
                
                % Core updraft with turbulence
                radius_variation = params.core_radius * ...
                    (1 + params.turbulence * sin(atan2(dy, dx) * 3) * ...
                    sin(height / params.noise_scale));
                core_updraft = params.strength * ...
                    (1 + params.turbulence * ...
                    (sin(dx / params.noise_scale) + ...
                     sin(dy / params.noise_scale) + ...
                     sin(height / params.noise_scale)) / 3) * ...
                    exp(-(drifted_r / radius_variation)^2);
                
                % Downdraft ring
                sink_profile = exp(-(r / params.sink_radius)^2) - ...
                               exp(-(r / params.core_radius)^2);
                sink_downdraft = params.strength * params.sink_strength * sink_profile;
                
                % Twisting component (optional turbulence)
                angle = atan2(dy, dx);
                Wx_rot = -dy/r * sin(height / params.noise_scale) * params.horiz_strength;
                Wy_rot = dx/r * sin(height / params.noise_scale) * params.horiz_strength;
                
                % Combine all components
                wind_field.Wx(j,i,k) = Wx_rot + params.turbulence * ...
                    sin(dx / params.noise_scale + height / params.noise_scale) + ...
                    params.constant_wind(1);
                wind_field.Wy(j,i,k) = Wy_rot + params.turbulence * ...
                    cos(dy / params.noise_scale + height / params.noise_scale) + ...
                    params.constant_wind(2);
                wind_field.U(j,i,k) = core_updraft + sink_downdraft + ...
                    params.constant_wind(3);
            end
        end
    end
end

% Add this after the other local function definitions:
function wind_field = create_simple_thermal_system(wind_field, params)
    % Set random seed for reproducible thermal positions
    rng(params.rng_seed);
    
    % Get grid dimensions
    x_size = max(wind_field.x) - min(wind_field.x);
    y_size = max(wind_field.y) - min(wind_field.y);
    
    % Generate random thermal positions
    x_pos = min(wind_field.x) + rand(params.num_thermals, 1) * x_size;
    y_pos = min(wind_field.y) + rand(params.num_thermals, 1) * y_size;
    
    % Generate random thermal parameters within specified ranges
    strengths = params.strength_range(1) + rand(params.num_thermals, 1) * ...
                (params.strength_range(2) - params.strength_range(1));
    radii = params.radius_range(1) + rand(params.num_thermals, 1) * ...
            (params.radius_range(2) - params.radius_range(1));
    
    % Initialize wind field
    wind_field.Wx = ones(size(wind_field.X)) * params.constant_wind(1);
    wind_field.Wy = ones(size(wind_field.Y)) * params.constant_wind(2);
    wind_field.U = ones(size(wind_field.U)) * params.constant_wind(3);
    
    % Add contribution from each thermal
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                % Sum contributions from all thermals
                for t = 1:params.num_thermals
                    r = sqrt((wind_field.X(j,i,k) - x_pos(t))^2 + ...
                            (wind_field.Y(j,i,k) - y_pos(t))^2);
                    % Add Gaussian thermal profile
                    wind_field.U(j,i,k) = wind_field.U(j,i,k) + ...
                        strengths(t) * exp(-(r/radii(t))^2);
                end
            end
        end
    end
end


function wind_field = create_drifting_thermal_system(wind_field, params)
    % Set random seed for reproducible thermal positions
    rng(params.rng_seed);
    
    % Get grid dimensions
    x_size = max(wind_field.x) - min(wind_field.x);
    y_size = max(wind_field.y) - min(wind_field.y);
    
    % Generate random thermal positions
    x_pos = min(wind_field.x) + rand(params.num_thermals, 1) * x_size;
    y_pos = min(wind_field.y) + rand(params.num_thermals, 1) * y_size;
    
    % Generate random thermal parameters within specified ranges
    strengths = params.strength_range(1) + rand(params.num_thermals, 1) * ...
                (params.strength_range(2) - params.strength_range(1));
    radii = params.radius_range(1) + rand(params.num_thermals, 1) * ...
            (params.radius_range(2) - params.radius_range(1));
    
    % Initialize wind field
    wind_field.Wx = ones(size(wind_field.X)) * params.constant_wind(1);
    wind_field.Wy = ones(size(wind_field.Y)) * params.constant_wind(2);
    wind_field.U = ones(size(wind_field.U)) * params.constant_wind(3);
    
    % Add contribution from each thermal
    for i = 1:length(wind_field.x)
        for j = 1:length(wind_field.y)
            for k = 1:length(wind_field.z)
                % Calculate height for drift calculation
                height = wind_field.Z(j,i,k);
                
                % Skip if above maximum thermal height
                if height > params.max_thermal_height
                    continue;  % No thermal contribution above max height
                end
                
                % Calculate height factor (1 at ground, 0 at max height)
                % This creates a smooth fade-out as thermals approach their max height
                height_factor = max(0, 1 - (height / params.max_thermal_height)^2);
                
                % Sum contributions from all thermals
                for t = 1:params.num_thermals
                    % Calculate drift based on height
                    drift_x = params.drift_direction(1) * height * params.drift_rate;
                    drift_y = params.drift_direction(2) * height * params.drift_rate;
                    
                    % Calculate distance from drifted thermal center
                    r = sqrt((wind_field.X(j,i,k) - (x_pos(t) + drift_x))^2 + ...
                            (wind_field.Y(j,i,k) - (y_pos(t) + drift_y))^2);
                    
                    % Add Gaussian thermal profile with height attenuation
                    wind_field.U(j,i,k) = wind_field.U(j,i,k) + ...
                        strengths(t) * exp(-(r/radii(t))^2) * height_factor;
                end
            end
        end
    end
end