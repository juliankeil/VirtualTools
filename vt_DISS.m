function out = vt_DISS(cfg,data1,data2)

% Computes Global Dissimilarity following the formula from Murray et al.,
% 2008.
% Requires two Timelock datasets as input.
% Attention! Needs to be Log-Transformed to test against 0!
%
% Use as:
% cfg.latency = latency for which to compute the GFP
% 
% out = vt_DISS(cfg,data1,data2);
%
% Ver 1.0.: 20.11.2015 Julian Keil
%% 1. Check input

t1 = nearest(data1.time,cfg.latency(1));
t2 = nearest(data1.time,cfg.latency(end));

timevec = t1:1:t2;
%% 2. Compute Dissimilarity

for t = 1:length(timevec)

    GFP1 = sqrt((1/length(data1.label)) .* sum(data1.avg(:,timevec(t)).^2));
    GFP2 = sqrt((1/length(data2.label)) .* sum(data2.avg(:,timevec(t)).^2));
    
    DISS(t) = sqrt((1/length(data1.label)) .* sum(((data1.avg(:,t)/GFP1)-(data2.avg(:,t)/GFP2)).^2));

end

%% 3. Create Output
out.label = {'diss'};
out.time = data1.time(timevec);
out.dimord = 'chan_time';
out.cfg = data1.cfg;
out.diss = DISS;
end