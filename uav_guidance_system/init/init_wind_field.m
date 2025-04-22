function wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing)
    % Creates a 3D grid structure for wind field data
    %
    % Inputs:
    %   x_limits = [x_min, x_max]
    %   y_limits = [y_min, y_max]
    %   z_limits = [z_min, z_max]
    %   grid_spacing = distance between grid points (m)
    %
    % Output:
    %   wind_field structure containing grid and wind data

    % Create grid coordinates
    x = x_limits(1):grid_spacing:x_limits(2);
    y = y_limits(1):grid_spacing:y_limits(2);
    z = z_limits(1):grid_spacing:z_limits(2);
    
    % Store grid information
    wind_field.x = x;
    wind_field.y = y;
    wind_field.z = z;
    wind_field.grid_spacing = grid_spacing;
    
    % Create 3D matrices for wind components
    [X, Y, Z] = meshgrid(x, y, z);
    wind_field.X = X;  % Grid X coordinates
    wind_field.Y = Y;  % Grid Y coordinates
    wind_field.Z = Z;  % Grid Z coordinates
    
    % Initialize wind component matrices
    wind_field.Wx = zeros(size(X));  % x-component of wind
    wind_field.Wy = zeros(size(Y));  % y-component of wind
    wind_field.U  = zeros(size(Z));  % z-component of wind (vertical)
end