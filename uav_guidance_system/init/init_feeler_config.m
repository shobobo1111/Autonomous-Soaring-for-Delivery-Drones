% init_feeler_config.m
function feeler_config = init_feeler_config()
% INIT_FEELER_CONFIG Initialize configuration for feeler line path extension
% Contains parameters controlling the behavior of look-ahead path evaluation

% Basic feeler control
feeler_config.enabled = true;         % Enable/disable feeler lines
feeler_config.length = 100;            % Length of feeler lines (meters)
feeler_config.num_points = 10;        % Number of evaluation points along feeler

% Reward weighting parameters
feeler_config.initial_weight = 1.0;   % Base weight for feeler rewards
feeler_config.decay_rate = 0.9;       % Decay rate for rewards along feeler

end