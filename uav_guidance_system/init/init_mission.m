% function mission = init_mission()
%     mission.time_to_deliver = 300; % time delivery should be made in in seconds - made larger for long dist opt
%     mission.alt_floor = 150; % min altitude during delivery, not used for take off and landing
%     mission.alt_floor_tol = 20;
% end


function mission = init_mission(varargin)
    % INIT_MISSION Initialize mission parameters with optional customization
    %
    % Usage:
    %   mission = init_mission()           % Return default mission params
    %   mission = init_mission(time_value) % Override time_to_deliver
    
    % Default mission parameters
    mission.time_to_deliver = 100;    % Time allocated for delivery (seconds)
    mission.alt_floor = 150;          % Minimum safe altitude (meters)
    mission.alt_floor_tol = 50;       % Tolerance for altitude floor warnings (meters)
    
    % Override time_to_deliver if provided
    if nargin > 0 && ~isempty(varargin{1})
        mission.time_to_deliver = varargin{1};
    end
end