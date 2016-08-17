function [out] = vt_entropy(cfg,indat)

% Compute the time-resolved entropy by sliding a 256-element wide window
% over the data.
% Window is centered around each time point, i.e. out(i) = entropy(in(i-77):in(i+78)) 
% 
% Result will be a n-256 element single trial estimation of entropy.
%
% Input
% cfg.parameter = Parameter to work on, default is trial
%
% (c) Julian Keil 24.06.2014
% Version 1: Inspired by Werkle-Berger et al. 2014

%%

if isfield(cfg,'parameter')
    parameter = cfg.parameter;
else
    parameter = 'trial';
end

out = indat; % Predefine Output

%% Loop Trials if present

if isfield(indat, 'trial')
    tmp = double(NaN(size(indat.(parameter){1}))); % Preallocate for Speed
    for t = 1:length(indat.(parameter)) % Loop Trials
        fprintf('Estimating Chaos at Trial %d\n',t)
        for c = 1:size(indat.(parameter){t},1) % Loop Channels
            for s = 128:size(indat.(parameter){t},2)-128 % Loop Time Points
                tmp(c,s) = entropy(double(indat.(parameter){t}(c,s-127:s+128)));
            end % s
        end % c
        out.trial{t} = tmp;
    end % t
    
else % Skip the trial part
    tmp = double(NaN(size(indat.(parameter)))); % Preallocate for Speed
     fprintf('Estimating Average Chaos %d\n')
        for c = 1:size(indat.(parameter),1) % Loop Channels
            for s = 128:size(indat.(parameter),2)-128 % Loop Time Points
                tmp(c,s) = entropy(double(indat.(parameter)(c,s-127:s+128)));
            end % s
        end % c
        out.(parameter) = tmp;
        out = rmfield(out,'var');
end
%%

