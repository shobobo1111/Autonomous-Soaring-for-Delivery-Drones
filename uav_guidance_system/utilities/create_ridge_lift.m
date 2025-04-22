
function wind_field = create_ridge_lift(terrain_data, wind_speed, wind_direction, wind_field)
    % CREATE_RIDGE_LIFT Update an existing 3D wind field with ridge lift effects
    %
    % Inputs:
    %   terrain_data: Structure containing terrain data and slope vectors
    %   wind_speed: Wind speed in m/s
    %   wind_direction: Wind direction in degrees (meteorological convention)
    %   wind_field: Existing wind field structure to modify (REQUIRED)
    %
    % Output:
    %   wind_field: Structure with 3D wind field components updated with ridge lift
    
    % Check if wind_field was provided
    if nargin < 4 || isempty(wind_field)
        error('Wind field must be provided. Initialize using init_wind_field() first.');
    end
    
    % Get dimensions of terrain and wind field
    [terrain_ny, terrain_nx] = size(terrain_data.terrain);
    [wind_ny, wind_nx, wind_nz] = size(wind_field.X);
    
    % Determine the overlapping dimensions to use
    ny = min(terrain_ny, wind_ny);
    nx = min(terrain_nx, wind_nx);
    nz = wind_nz;
    
    if ny ~= terrain_ny || nx ~= terrain_nx || ny ~= wind_ny || nx ~= wind_nx
        warning(['Dimension mismatch between terrain [%d, %d] and wind field [%d, %d]. ', ...
                'Using smaller dimensions [%d, %d] for both.'], ...
                terrain_ny, terrain_nx, wind_ny, wind_nx, ny, nx);
    end
    
    fprintf('Processing area: [%d, %d]\n', ny, nx);
    fprintf('Terrain elevation range: %.2f to %.2f m\n', min(min(terrain_data.terrain(1:ny,1:nx))), max(max(terrain_data.terrain(1:ny,1:nx))));
    
    % Convert wind direction to vector components (meteorological to mathematical)
    wind_dir_rad = (270 - wind_direction) * pi/180;
    wind_unit_vector = [cos(wind_dir_rad), sin(wind_dir_rad)];  % Unit vector in wind direction
    
    % Height-based decay parameters
    decay_height = 150;  % Height where effect reduces to ~37% (meters)
    min_effect = 0.0;    % Minimum terrain effect at high altitudes

    % % Initialise wind field to ambient conditions (important for below-terrain points)
    % for k = 1:nz
    %     wind_field.Wx(:,:,k) = wind_speed * cos(wind_dir_rad);
    %     wind_field.Wy(:,:,k) = wind_speed * sin(wind_dir_rad);
    %     wind_field.U(:,:,k) = 0; % No vertical component by default
    % end
    
    % Process each grid point in the overlapping region
    for i = 1:ny
        for j = 1:nx
            % Get terrain height at this location
            terrain_height = terrain_data.terrain(i,j);
            
            % Get terrain slope vectors
            slope_x = terrain_data.slope_dx(i,j);
            slope_y = terrain_data.slope_dy(i,j);
            slope_z = terrain_data.slope_dz(i,j);  % Vertical component
            
            % Project to horizontal plane for angle calculation
            slope_horiz = [slope_x, slope_y];
            slope_horiz_mag = norm(slope_horiz);
            
            % If slope is too flat, skip detailed calculations
            if slope_horiz_mag < 0.001
                continue; % Already initialized to ambient wind
            end
            
            % Normalize horizontal component
            slope_horiz_unit = slope_horiz / slope_horiz_mag;
            
            % Calculate dot product between wind and upslope direction
            dot_product = dot(wind_unit_vector, slope_horiz_unit);
            
            % Calculate absolute value of vertical slope component (steepness)
            steepness = abs(slope_z);
            
            % Combined factor: alignment * steepness
            % Only create strong upward deflection when wind blows INTO the slope
            if dot_product > 0  % Wind blowing into slope
                deflection_factor = dot_product * steepness;
            else  % Wind blowing away from slope - minimal lift or slight down
                deflection_factor = dot_product * steepness * 0.1;  % Reduced effect
            end
            
            % Process each height level with decay
            for k = 1:nz
                % Get absolute height of this grid point
                absolute_height = wind_field.z(k);
                fprintf('%.f',absolute_height);
                
                % Only process points above the terrain
                if absolute_height >= terrain_height
                    % Calculate height above ground level
                    height_agl = absolute_height - terrain_height;
                    
                    % Calculate decay factor based on height above ground
                    decay_factor = max(min_effect, ...
                                    (1 - min_effect) * exp(-height_agl/decay_height) + min_effect);
                    
                    % Apply horizontal wind components (unchanged)
                    wind_field.Wx(i,j,k) = wind_speed * cos(wind_dir_rad);
                    wind_field.Wy(i,j,k) = wind_speed * sin(wind_dir_rad);
                    
                    % Calculate vertical component based on slope interaction
                    vertical_component = deflection_factor * wind_speed * decay_factor;
                    wind_field.U(i,j,k) = vertical_component;
                else
                    % For points below terrain, set to zero (or you could implement
                    % a "no-flow" boundary condition by setting all components to zero)
                    wind_field.Wx(i,j,k) = 0;
                    wind_field.Wy(i,j,k) = 0;
                    wind_field.U(i,j,k) = 0;
                end
            end
        end
    end
    
    fprintf('Wind field updated with ridge lift effects accounting for terrain height\n');
end



%     % Process each grid point in the overlapping region
%     for i = 1:ny
%         for j = 1:nx
%             % Get terrain slope vectors from the new data format
%             slope_x = terrain_data.slope_dx(i,j);
%             slope_y = terrain_data.slope_dy(i,j);
%             slope_z = terrain_data.slope_dz(i,j);  % Vertical component
% 
%             % Project to horizontal plane for angle calculation
%             slope_horiz = [slope_x, slope_y];
%             slope_horiz_mag = norm(slope_horiz);
% 
%             % If slope is too flat, just use ambient wind
%             if slope_horiz_mag < 0.001
%                 for k = 1:nz
%                     wind_field.Wx(i,j,k) = wind_speed * cos(wind_dir_rad);
%                     wind_field.Wy(i,j,k) = wind_speed * sin(wind_dir_rad);
%                     wind_field.U(i,j,k) = 0;
%                 end
%                 continue;
%             end
% 
%             % Normalize horizontal component
%             slope_horiz_unit = slope_horiz / slope_horiz_mag;
% 
%             % Calculate dot product between wind and upslope direction
%             % This gives cosine of angle between vectors
%             dot_product = dot(wind_unit_vector, slope_horiz_unit);
% 
%             % Calculate absolute value of vertical slope component
%             % This represents steepness
%             steepness = abs(slope_z);
% 
%             % Combined factor: alignment * steepness
%             % Only create strong upward deflection when wind blows INTO the slope
%             if dot_product > 0  % Wind blowing into slope
%                 deflection_factor = dot_product * steepness;
%             else  % Wind blowing away from slope - minimal lift or slight down
%                 deflection_factor = dot_product * steepness * 0.1;  % Reduced effect when flowing away
%             end
% 
%             % Process each height level with decay
%             for k = 1:nz
%                 height_agl = wind_field.z(k);
% 
%                 % Calculate decay factor based on height
%                 decay_factor = max(min_effect, ...
%                                 (1 - min_effect) * exp(-height_agl/decay_height) + min_effect);
% 
%                 % Apply horizontal wind components (unchanged)
%                 wind_field.Wx(i,j,k) = wind_speed * cos(wind_dir_rad);
%                 wind_field.Wy(i,j,k) = wind_speed * sin(wind_dir_rad);
% 
%                 % Calculate vertical component based on slope interaction
%                 % Stronger effect when wind aligns with upslope direction
%                 vertical_component = deflection_factor * wind_speed * decay_factor;
%                 wind_field.U(i,j,k) = vertical_component;
%             end
%         end
%     end
% 
%     fprintf('Wind field updated with ridge lift effects\n');
% end




% function wind_field = create_ridge_lift(terrain_data, wind_speed, wind_direction)
%     % CREATE_RIDGE_LIFT Generate a 3D wind field with ridge lift effects
%     %
%     % Inputs:
%     %   terrain_data: Structure containing terrain data and slope vectors
%     %   wind_speed: Wind speed in m/s
%     %   wind_direction: Wind direction in degrees (meteorological convention)
%     %
%     % Output:
%     %   wind_field: Structure with 3D wind field components
% 
%     % Print debug information about terrain data
%     fprintf('Terrain dimensions: [%d, %d]\n', size(terrain_data.terrain));
%     fprintf('Z_AGL points: %d\n', length(terrain_data.z_agl));
%     fprintf('Terrain elevation range: %.2f to %.2f m\n', min(min(terrain_data.terrain)), max(max(terrain_data.terrain)));
% 
%     % Create coordinate grid
%     % Convert lat/lon to meters for easier visualization
%     lat_center = mean(terrain_data.y);
%     lon_scale = 111000 * cos(lat_center * pi/180);  % Approximate meters per degree longitude
%     lat_scale = 111000;  % Approximate meters per degree latitude
% 
%     % Create coordinate vectors in meters
%     x_meters = (terrain_data.x - mean(terrain_data.x)) * lon_scale;
%     y_meters = (terrain_data.y - mean(terrain_data.y)) * lat_scale;
% 
%     % Create 3D grid
%     [X, Y, Z] = ndgrid(y_meters, x_meters, terrain_data.z_agl);
% 
%     % Store grid coordinates for wind field
%     wind_field.x = x_meters;
%     wind_field.y = y_meters;
%     wind_field.z = terrain_data.z_agl;
%     wind_field.X = X;
%     wind_field.Y = Y;
% 
%     % Get terrain heights and mesh with Z coordinates
%     terrain_heights = terrain_data.terrain;
%     [ny, nx] = size(terrain_heights);
%     nz = length(terrain_data.z_agl);
% 
%     % Create terrain surface
%     Z_terrain = zeros(size(X));
%     for k = 1:nz
%         Z_terrain(:,:,k) = repmat(terrain_heights, [1, 1]) + terrain_data.z_agl(k);
%     end
%     wind_field.Z = Z_terrain;
% 
%     % Convert wind direction to vector components (meteorological to mathematical)
%     wind_dir_rad = (270 - wind_direction) * pi/180;
%     wind_unit_vector = [cos(wind_dir_rad), sin(wind_dir_rad)];  % Unit vector in wind direction
% 
%     % Initialize wind components
%     wind_field.Wx = zeros(size(X));
%     wind_field.Wy = zeros(size(Y));
%     wind_field.U = zeros(size(Z));
% 
%     % Height-based decay parameters
%     decay_height = 150;  % Height where effect reduces to ~37% (meters)
%     min_effect = 0.0;    % Minimum terrain effect at high altitudes
% 
%     % Process each grid point
%     for i = 1:ny
%         for j = 1:nx
%             % Get terrain slope vectors from the new data format
%             slope_x = terrain_data.slope_dx(i,j);
%             slope_y = terrain_data.slope_dy(i,j);
%             slope_z = terrain_data.slope_dz(i,j);  
% 
%             % Project to horizontal plane for angle calculation
%             slope_horiz = [slope_x, slope_y];
%             slope_horiz_mag = norm(slope_horiz);
% 
%             % If slope is too flat, just use ambient wind
%             if slope_horiz_mag < 0.001
%                 for k = 1:nz
%                     wind_field.Wx(i,j,k) = wind_speed * cos(wind_dir_rad);
%                     wind_field.Wy(i,j,k) = wind_speed * sin(wind_dir_rad);
%                     wind_field.U(i,j,k) = 0;
%                 end
%                 continue;
%             end
% 
%             % Normalize horizontal component
%             slope_horiz_unit = slope_horiz / slope_horiz_mag;
% 
%             % Calculate dot product between wind and upslope direction
%             % This gives cosine of angle between vectors
%             dot_product = dot(wind_unit_vector, slope_horiz_unit);
% 
%             % Calculate absolute value of vertical slope component
%             % This represents steepness
%             steepness = abs(slope_z);
% 
%             % Combined factor: alignment * steepness
%             deflection_factor = dot_product * steepness;
% 
%             % Process each height level with decay
%             for k = 1:nz
%                 height_agl = terrain_data.z_agl(k);
% 
%                 % Calculate decay factor based on height
%                 decay_factor = max(min_effect, ...
%                                 (1 - min_effect) * exp(-height_agl/decay_height) + min_effect);
% 
%                 % Apply horizontal wind components (unchanged)
%                 wind_field.Wx(i,j,k) = wind_speed * cos(wind_dir_rad);
%                 wind_field.Wy(i,j,k) = wind_speed * sin(wind_dir_rad);
% 
%                 % Calculate vertical component based on slope interaction
%                 % Stronger effect when wind aligns with upslope direction
%                 vertical_component = deflection_factor * wind_speed * decay_factor;
%                 wind_field.U(i,j,k) = vertical_component;
%             end
%         end
%     end
% 
%     fprintf('Wind field generated: [%d, %d, %d]\n', size(wind_field.Wx));
% end

% function wind_field = create_ridge_lift(terrain_data, wind_speed, wind_direction)
%     % Debug prints
%     fprintf('Terrain data dimensions: [%d, %d]\n', size(terrain_data.terrain));
%     fprintf('Original z_agl points: %d\n', length(terrain_data.z_agl));
%     fprintf('Size of dx: [%s]\n', num2str(size(terrain_data.dx)));
%     fprintf('Size of dy: [%s]\n', num2str(size(terrain_data.dy)));
%     fprintf('Size of dz: [%s]\n', num2str(size(terrain_data.dz)));
%     fprintf('Range of terrain: %.2f to %.2f\n', min(min(terrain_data.terrain)), max(max(terrain_data.terrain)));
%     fprintf('Range of z_agl: %.2f to %.2f\n', min(terrain_data.z_agl), max(terrain_data.z_agl));
%     fprintf('X range: %.2f to %.2f\n', min(terrain_data.x), max(terrain_data.x));
%     fprintf('Y range: %.2f to %.2f\n', min(terrain_data.y), max(terrain_data.y));
% 
%     % Convert lat/lon to meters
%     lat_center = mean(terrain_data.y);
%     lon_scale = 111000 * cos(lat_center * pi/180);
%     lat_scale = 111000;
%     x_meters = (terrain_data.x - mean(terrain_data.x)) * lon_scale;
%     y_meters = (terrain_data.y - mean(terrain_data.y)) * lat_scale;
% 
%     % Create exact grid matching terrain data dimensions
%     x_coarse = linspace(min(x_meters), max(x_meters), length(terrain_data.x));
%     y_coarse = linspace(min(y_meters), max(y_meters), length(terrain_data.y));
% 
%     % Use existing height levels from terrain_data
%     agl_levels = terrain_data.z_agl;
% 
%     % % Initialise base meshgrid
%     % [X, Y] = meshgrid(x_coarse, y_coarse);
% 
%     % Create ndgrid instead of meshgrid
%     [X, Y, Z] = ndgrid(y_coarse, x_coarse, agl_levels);  % Note swapped x,y order
% 
% 
%     % Create 3D arrays - replicate X and Y for each height level
%     [ny, nx, nz] = size(X);
%     fprintf('Grid dimensions: [%d, %d, %d]\n', ny, nx, nz);
%     % nz = length(agl_levels);
% 
%     % % Create full 3D coordinate grids
%     % wind_field.X = repmat(X, [1, 1, nz]);
%     % wind_field.Y = repmat(Y, [1, 1, nz]);
%     % wind_field.Z = zeros(ny, nx, nz);
% 
%     % Store the coordinate matrices
%     wind_field.X = X;
%     wind_field.Y = Y;
%     wind_field.Z = Z;
% 
%     % Fill Z coordinates accounting for terrain
%     terrain_heights = terrain_data.terrain;
%     for k = 1:nz
%         wind_field.Z(:,:,k) = terrain_heights + agl_levels(k);
%     end
% 
%     % Store coordinate vectors for interpolation 
%     wind_field.x = x_coarse;
%     wind_field.y = y_coarse;
%     wind_field.z = agl_levels;
% 
%     fprintf('Wind field size: [%d, %d, %d]\n', ny, nx, nz);
% 
%     % Convert wind direction to vector components
%     wind_dir_rad = (270 - wind_direction) * pi/180;
%     wind_horizontal = [cos(wind_dir_rad), sin(wind_dir_rad), 0];  % Unit vector in wind direction
%     wind_vector = wind_horizontal * wind_speed;
% 
%     % Initialise the wind components
%     wind_field.Wx = zeros(size(X));
%     wind_field.Wy = zeros(size(Y));
%     wind_field.U = zeros(size(Z));
% 
%     % Height-based decay parameters
%     decay_height = 750;  % Height over which effect decays to negligible (meters) - TO TUNE
%     min_effect = 0.0;    % Minimum terrain effect factor - TO TUNE
% 
%     % For each point above terrain
%     for i = 1:ny
%         for j = 1:nx
%             terrain_height = terrain_heights(i,j);
% 
%             for k = 1:nz
%                 if wind_field.Z(i,j,k) > terrain_height
%                     % Calculate height above terrain for this point
%                     height_agl = wind_field.Z(i,j,k) - terrain_height;
% 
%                     % Calculate decay factor (1 at surface, decreasing with height)
%                     decay_factor = max(min_effect, ...
%                         (1 - min_effect) * exp(-height_agl/decay_height) + min_effect);
% 
%                     % Get slope vector at this point (normalized)
%                     slope_vector = [terrain_data.dx(i,j,k), ...
%                                   terrain_data.dy(i,j,k), ...
%                                   terrain_data.dz(i,j,k)];
% 
%                     % Project and mirror horizontal component
%                     horiz_component = [slope_vector(1), slope_vector(2), 0];
%                     dot_product = dot(horiz_component, wind_horizontal);
%                     mirrored_horiz = 2 * dot_product * wind_horizontal - horiz_component;
% 
%                     % Combine with vertical component
%                     terrain_effect = [mirrored_horiz(1), mirrored_horiz(2), slope_vector(3)];
% 
%                     % Blend between terrain effect and free-stream wind based on decay
%                     resultant = terrain_effect * decay_factor * wind_speed + ...
%                               wind_vector * (1 - decay_factor);
% 
%                     % Store components
%                     wind_field.Wx(i,j,k) = resultant(1);
%                     wind_field.Wy(i,j,k) = resultant(2);
%                     wind_field.U(i,j,k) = resultant(3);
%                 end
%             end
%         end
%     end
% end