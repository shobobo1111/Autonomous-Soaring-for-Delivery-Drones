% First let's examine the contents of the NC file
filename = 'surface_sensible_heat_flux.nc';  % Replace with your .nc filename

% Read main variables using the actual names from your file
lat = ncread(filename, 'latitude');
lon = ncread(filename, 'longitude');
temp_2m = ncread(filename, 't2m');         % 2m temperature
sensible_heat = ncread(filename, 'slhf');  % Surface latent heat flux

% Let's examine the data structure
fprintf('Data dimensions:\n');
fprintf('Latitude points: %d\n', length(lat));
fprintf('Longitude points: %d\n', length(lon));
fprintf('Temperature data size: [%s]\n', num2str(size(temp_2m)));

% Optional: Plot a sample temperature map to verify the data
figure;
pcolor(lon, lat, temp_2m(:,:,1)'); % Transpose for correct orientation
colorbar;
title('2m Temperature Map');
xlabel('Longitude');
ylabel('Latitude');