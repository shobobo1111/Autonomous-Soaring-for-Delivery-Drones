% plot_wind_field.m
function plot_wind_field(wind_field)
    % Create 3D quiver plot of entire wind field
    figure
    
    % Create 3D quiver plot (arrows showing wind vectors)
    quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
           wind_field.Wx, wind_field.Wy, wind_field.U, ...
           'b', 'LineWidth', 1);  % Made arrows thicker
    
    % Format plot
    xlabel('X (m)', 'FontSize', 12)  % Added font size
    ylabel('Y (m)', 'FontSize', 12)
    xlabel('Z (m)', 'FontSize', 12)
    title('3D Wind Field', 'FontSize', 14)  % Made title bigger
    grid on
    view(3)  % 3D view
    axis equal
    
    % Add colorbar to show wind strength
    wind_magnitude = sqrt(wind_field.Wx.^2 + wind_field.Wy.^2 + wind_field.U.^2);
    c = colorbar;
    ylabel(c, 'Wind Speed (m/s)')  % Added colorbar label
    colormap jet
end