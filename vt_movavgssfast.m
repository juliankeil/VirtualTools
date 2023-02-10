function [datamovavg,onsettimes]=vt_movavgssfast(cfg,varargin)
% Takes single trial data, computes the trial average and then 
% does the moving average of the number of cycles of the 
% specified modulation frequency
%
% cfg.trials = which trials to work on
% cfg.channel = which channels to work on
% cfg.modfreq = modulation frequency to extract
% cfg.cycnum = number of cycles to base the moving average on
% cfg.discard = how many cycles to discard, i.e. start point
% cfg.window = in what time window should the moving average be computed?
% cfg.baseline = [start stop] for baseline correction
%
%Nathan & Winny 2005
%Nathan 2006: deleted ampcrit and added preallocation for speeding up
% Julian 2023
%% 0. set the cfgs
modfreq = cfg.modfreq;
cycnum = cfg.cycnum;
discardcycles = cfg.discard;
samprat = varargin{1}.fsample;
wnd = cfg.window;

%% 1. Compute the ERP
icfg = [];
icfg.trials = cfg.trials;
icfg.channel = cfg.channel;

tmp_erp = ft_timelockanalysis(icfg,varargin{1});

% 1.1. Baseline correction
icfg = [];
icfg.baseline = cfg.baseline;

tmp_erp = ft_timelockbaseline(icfg,tmp_erp);

%% 2. compute the basic properties for the moving average
x = samprat/1000; % how many points per ms
startms = discardcycles*1000/modfreq; %a when to start the moving average
intervallms = 1000/modfreq; % one cycle in ms
intervallspt = intervallms*x; % one cylce in sampling points

startspt = floor(startms*x);
windowpoints = floor(cycnum*intervallspt);
maxend = startspt+ceil(wnd*samprat)-windowpoints;

% Prepare the output
datamovavg = tmp_erp;
datamovavg.time = tmp_erp.time(startspt:startspt+windowpoints-1);
datamovavg = rmfield(datamovavg,{'dof'});

% 2.1 first go to preallocate ressources

j=1;

startms2 = startms;
startspt2 = startspt;

while startspt2 < maxend
    startms2 = startms2+intervallms;
    startspt2 = floor(startms2*x);
    j = j+1;
end

onsettimes = zeros(j,2);
datatimecourse = zeros(size(tmp_erp.avg,1),windowpoints,j);

% 2.2 now real thing
i=1;

while startspt < maxend
    onsettimes(i,:) = [startms,startspt];

    datatimecourse(:,:,i) = tmp_erp.avg(:,startspt:startspt+windowpoints-1);
       
    startms = startms+intervallms;
    startspt = floor(startms*x);
    i = 1+i;
end

% Collect the output
datamovavg.avg = mean(datatimecourse,3);

