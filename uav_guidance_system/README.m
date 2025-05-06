% This project aims to find the optimal route for a fixed wing delivery
% drone given a wind feild

% To run: test_mutli_step_path.m

% /
% │
% ├── init/
% │   ├── init_aircraft.m       
% │   │   ├── Mass
% │   │   ├── Wing Dimensions
% │   │   ├── Drag coefficients
% │   │   ├── Lift coefficients
% │   │   ├── Speed limits
% │   │   └── Bank angle limits
% │   │
% │   ├── init_wind_field.m     
% │   │   ├── Creates grid coordinates
% │   │   ├── Creates 3D matrices for wind components
% │   │   └── Initialises wind components Wx, Wy, U 
% │   │
% │   ├── init_feeler_config.m     
% │   │   ├── Defines if enabled
% │   │   ├── Defines feeler length, No. eval points 
% │   │   └── Defines weight and decay for reward
% │   │
% │   └── init_waypoint.m 
% │       ├── Contains navigation time pressure parameters
% │       └── Defines and returns waypoint location
% │
% ├── core/
% │   ├── generate_paths.m      
% │   │   ├── Takes in current_state, aircraft, wind_field, elapsed_time
% │   │   ├── Gets init feeler_config and waypoint
% │   │   ├── Sets parameters:  bank, pitch, thrust, predict time, no.steps
% │   │   └── Generates paths, calculates rewards, sorts paths
% │   │
% │   ├── calculate_reward.m  
% │   │   ├── Potential energy (R_P)
% │   │   ├── Kinetic energy (R_K) 
% │   │   ├── Updraft reward (R_U)
% │   │   ├── Navigation reward (R_N) 
% │   │   └── Returns total reward calculation
% │   │
% │   ├── flight_dynamics.m
% │   │   ├── 3-DOF equations
% │   │   ├── State propagation
% │   │   └── Energy calculations
% │   │
% │   ├── generate_multi_step_path.m
% │   │   ├── takes in the inital position
% │   │   ├── calls the generate paths function
% │   │   │   └── returns list of paths from highest reward function 
% │   │   └── picks the path with highest reward and moves to next step
% │   │
% │   └── get_wind_at_position.m        
% │       ├── Interpolates wind components at any x,y,z position
% │       ├── Takes in (x, y, z, wind_field)
% │       └── Outputs: wind: Structure with Wx, Wy, U at requested position
% │       
% │
% ├── utilities/  
% │   ├── plot_wind_feild.m        
% │   │   ├── takes in wind_feild
% │   │   ├── visual plots wind feild vectors
% │   │   └── called in test_wind_feild.m
% │   │
% │   ├── plot_path_flower.m  
% │   │   ├── takes in paths, wind_field, step_number, elapsed_time
% │   │   ├── Outputs sert of possible paths visually
% │   │   └── called at each step within test_mutli_step_path.m
% │   │
% │   └──  create_thermals.m  
% │       ├── Takes in wind_feild and updraft type
% │       ├── Contains a number of different updrafts
% │       └── Returns requested updraft model
% │   
% ├── testing/
% │   ├── test_multi_step_path.m <------------------------------------------- RUN THIS FILE TO START
% │   │   ├── Requires user to pick updraft type and other parameters
% │   │   ├── Calls and tests the generate_multi_step_path.m
% │   │   └── Plots steps using plot_path_flower.m, plots full path in file
% │   │
% │   ├── test_path_generation
% │   │   ├── User sets basic parameters
% │   │   ├── Calls and tests the generate_paths.m
% │   │   └── Plots top down view and path flower in file
% │   │
% │   ├── test_flight_dynamics_2.m  
% │   │   ├── Test straight flight
% │   │   ├── Test turning flight
% │   │   ├── Test wind effects
% │   │   ├── Test energy conservation
% │   │   └── Visualise results
% │   │
% │   └── test_wind_field.m
% │       └── calls the init_wind_feild.m to test it
% │
% └── startup.m  
%     ├──projectPath = fileparts(mfilename('fullpath'));
%     └──addpath(genpath(projectPath));

