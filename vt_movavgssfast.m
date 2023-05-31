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
%% 0. set the cfgs and defaults
if isfield(cfg, 'modfreq')
    modfreq = cfg.modfreq;
else
    fprintf('Please set the target frequency \n')
    return
end

if isfield(cfg, 'cycnum')
    cycnum = cfg.cycnum;
else
    cycnum = 5;
    fprintf('No cycle number set, using the default of 5 cycles \n')
end

if isfield(cfg, 'discard')
    discardcycles = cfg.discard;
else
    discardcycles = 1;
    fprintf('No number of to-be-discarded cycles set, using the default of 1 cycle \n')
end

if isfield(cfg, 'window')
    wnd = cfg.window;
else
    fprintf('Please set the analysis window \n')
    return
end

if ~isfield(cfg, 'trials')
    cfg.trials = 'all';
    fprintf('Using all available trials \n')
end

if ~isfield(cfg, 'channel')
    cfg.channel = 'all';
    fprintf('Using all available channels \n')
end

if ~isfield(cfg, 'baseline')
    cfg.baseline = [varargin{1}.time{1}(1) varargin{1}.time{1}(end)];
    fprintf('Using the entire window for baseline correction \n')
    return
end

%% 1. Downsample to maximum multiples of sample rate
samprat = varargin{1}.fsample;
icfg = [];
icfg.resamplefs = modfreq*(floor(samprat/modfreq));

tmp_dat = ft_resampledata(icfg,varargin{1});
samprat = tmp_dat.fsample; % Update Sample rate

%% 2. Compute the average across trials
icfg = [];
icfg.trials = cfg.trials;
icfg.channel = cfg.channel;

tmp_erp = ft_timelockanalysis(icfg,tmp_dat);

% 2.1. Baseline correction
icfg = [];
icfg.baseline = cfg.baseline;

tmp_erp = ft_timelockbaseline(icfg,tmp_erp);

%% 3. compute the basic properties for the moving average
mspts = samprat/1000; % how many points per ms
startms = discardcycles*1000/modfreq; %a when to start the moving average
intervallms = 1000/modfreq; % one cycle in ms
intervallspt = intervallms*mspts; % one cylce in sampling points

startspt = floor(startms*mspts);
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
    startspt2 = floor(startms2*mspts);
    j = j+1;
end

onsettimes = NaN(j,2);
datatimecourse = NaN(size(tmp_erp.avg,1),windowpoints,j);

% 2.2 now real thing
i=1;

while startspt < maxend
    onsettimes(i,:) = [startms,startspt];

    datatimecourse(:,:,i) = tmp_erp.avg(:,startspt:startspt+windowpoints-1);
       
    startms = startms+intervallms;
    startspt = floor(startms*mspts);
    i = 1+i;
end

% Collect the output
datamovavg.avg = mean(datatimecourse,3,'omitnan');
datamovavg.fsample = samprat;

