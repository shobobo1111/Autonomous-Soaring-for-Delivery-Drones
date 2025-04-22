function waypoint = init_waypoint()
    % % Position
    waypoint.x = -100;
    waypoint.y = 100;
    waypoint.z = 600;

    waypoint.arrival_tol = 50; % within this distance (meters) the program considers itself arrived

    waypoint.loiter_within_r = 400; % this changes the waypoint to a loitering area, defined by radius in meters

    % % Set a distant waypoint to establish direction for range testing
    % waypoint.x = 0;  % 5 km away in x-direction
    % waypoint.y = 5000;  % 5 km away in y-direction
    % waypoint.z = 450;   % Same altitude as initial
    
end