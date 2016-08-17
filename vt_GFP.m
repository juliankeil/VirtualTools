function out = vt_GFP(cfg,data)

% Computes Global Field Power following the formula from Murray et al.,
% 2008.
% Requires Timelock data as input.
%
% Use as:
% cfg.latency = latency for which to compute the GFP
% 
% out = vt_GFP(cfg,data);
%
% Ver 1.0.: 20.11.2015 Julian Keil
%% 1. Check input

t1 = nearest(data.time,cfg.latency(1));
t2 = nearest(data.time,cfg.latency(end));

timevec = t1:1:t2;
%% 2. Compute GFP

for t = 1:length(timevec)

    GFP(t) = sqrt((1/length(data.label)) .* sum(data.avg(:,timevec(t)).^2));
    
end

%% 3. Create Output
out.label = {'gfp'};
out.time = data.time(timevec);
out.dimord = 'chan_time';
out.cfg = data.cfg;
out.gfp = GFP;
end