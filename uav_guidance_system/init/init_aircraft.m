function aircraft = init_aircraft()
    % AIRCRAFT_CONFIG Returns configuration for 10kg fixed wing UAV
    % Physical Properties
    aircraft.mass = 9;          % kg
    aircraft.wingspan = 3;       % m
    aircraft.wing_area = 1;      % m^2
    aircraft.chord = aircraft.wing_area/aircraft.wingspan;  % m, assuming rectangular wing
    
    % Aerodynamic Coefficients
    aircraft.CL0 = 0.3;          % Base lift coefficient (UPDATE LATER with real data)
    aircraft.CL_alpha = 5.7;
    aircraft.CD0 = 0.015;        % Zero-lift drag coefficient (UPDATE LATER)
    aircraft.k = 0.015;          % Induced drag factor (UPDATE LATER)
    
    % Control Surface Properties
    aircraft.elevator_max = 30 * (pi/180);  % max elevator deflection (rad)
    aircraft.elevator_effect = 2.0;         % elevator effectiveness (UPDATE LATER)
    
    % Inertial Properties
    aircraft.I_yy = 0.5;         % Pitch moment of inertia (UPDATE LATER)
    
    % Trim Conditions (for reference)
    aircraft.cruise_speed = 15;   % m/s
    aircraft.min_speed = 10;      % m/s
    aircraft.max_speed = 30;      % m/s
    aircraft.max_climb_angle = 15 * (pi/180);  % 15 degrees max climb angle
end