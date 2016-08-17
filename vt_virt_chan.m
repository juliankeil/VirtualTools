function [sp] = vt_virt_chan(cfg,alldata,indata);
% Use as: projected_data = vt_virt_chan(cfg,alldata,indata);
%
% Builds a spatial filter for locations of interest from "alldata" and
% projects "indata" through filter.
% Filter can be time-frequency-optimized (DICS) or broadband (LCMV)
%
% cfg.location= location in Head Corrdinates e.g. [52 -53 4; 12 -37 24];
% OR
% cfg.lf = precomputed leadfield
% cfg.vol= Head Volume Model;
% cfg.elec= Electrode Structure;
%
% cfg.beamformer= Method for the spatial filter'dics' or 'lcmv'
% cfg.lambda= specify regularization parameter, e.g. '5%', 0.05;
%
% In case of DICS:
%   cfg.foi= Frequency of Interest for Filter;
%   cfg.toi= Time-Window of Interest for Filter
%   cfg.tapsmofrq= taper smoothing freuency in Hz;
%   cfg.taper = Taper Window, default is "Hanning"
%   cfg.csdmethod = either "fft" or "tfr"
%   if you select "tfr" you need to specify
%       cfg.t_ftimwin = Time-Frequency Window length .4;
%       cfg.foilim = Frequency Range
%       cfg.latency = Time of Interest to center Time-Frequency Window
% In case of LCMV:
%   cfg.covwin = Window for the Computation of Covariance, taken of alldat -> The longer the
%   better

%
% (c) Julian Keil 2013
% Version 1: 
% Version 2: 03.09.2013
% Version 2.1: 16.09.2013 edited by Martin Krebber, minor changes in 'Use
% as:' and 4.1 DICS (method='mtmfft', latency no longer required)
% Version 2.2: 12.11.2013: Added MEG-compatibility and automatic
% electrode-stucture generation
% Version 2.3: 15.11.2013: Added cutting of trials to optimize DICS-filter
% Version 2.4: 18.11.2013: Added TFR and Taper Options
% Version 2.5: 16.01.2014: Added Leadfield Normalization Option. Fixed
% Label generation for the "location"-option

%% 1. Set Basics

vol = cfg.vol;
elec = cfg.elec;
label = alldata.label;
    elec.label=label;
    
if isfield(cfg,'location')
    location = cfg.location;
    tmplf = [];
elseif isfield(cfg,'lf')
    tmplf = cfg.lf;
    location = [];
else
    disp('please provide either leadfield or locations')
    return
end

if isfield(cfg,'lambda')
    lambda = cfg.lambda;
else
    lambda = 0;
end

method = cfg.beamformer;

if strcmpi(method,'dics')
    if isfield(cfg,'foi')
        foi = cfg.foi;
    elseif isfield(cfg,'foilim')
        foi = mean(foilim);
    else
        disp('please state either FOI or FOILIM')
        return
    end
    
    if isfield(cfg,'latency')
        latency = cfg.latency;
    else
        latency = alldata.time{1}(ceil(length(alldata.time{1})/2));
    end
    
    if isfield(cfg,'csdmethod')
        csdmethod = cfg.csdmethod;
    else
        csdmethod = 'fft';
    end
    
    if strcmpi(csdmethod,'tfr')
        t_ftimwin = cfg.t_ftimwin;
    end
    
    tapsmofrq = cfg.tapsmofrq;
    
    if isfield(cfg,'taper')
        taper = cfg.taper;
    else
        taper = 'hanning';
    end
    
elseif strcmpi(method,'lcmv');
    covwin = cfg.covwin;
else
    disp('please state beamformer method')
    return
end

% Leadfield Normalize Options
if isfield(cfg,'normalize')
    normalize = cfg.normalize;
else
    normalize = 'no';
end

%% 2. Build new Leadfield specific for our location
if isempty(tmplf);
    disp('Building a new Leadfield')
    
    vol.unit='mm';
    elec.unit='mm';

    tmpcfg = [];
    %cfg.channel = alldata.label; % all EEG Channels
    tmpcfg.elec = elec;
    tmpcfg.vol = vol;
    
    % Set Grid Points
    tmpcfg.grid.pos=location; % Coordinates -> Look up in the source plot
    tmpcfg.grid.inside=1:size(location,1);
    
    % Leadfield Normalization Option
    if strcmpi(normalize,'yes')
        tmpcfg.normalize = 'yes';
    else
        tmpcfg.normalize = 'no';
    end
    
    lf2=ft_prepare_leadfield(tmpcfg);
    lf2.unit='mm';
    lf2.cfg.elec.label=label;
    
else
    lf2 = tmplf;
    lf2.unit = 'mm';
    vol.unit = 'mm';
    elec.unit= 'mm';
end

%% 3. Make a new Time Vector for all input data sets
timevec =  alldata.time{1}(1):(1/alldata.fsample):alldata.time{1}(end);
intime = indata.time{1};

for t=1:length(alldata.trial);
    alldata.time{t} = timevec;
end


for t=1:length(indata.trial);
    indata.time{t} = timevec;
end

%% 4. Switch between DICS and LCMV

switch method
    case{'dics'}
        disp('Building the DICS Filter')
	
        if isfield(cfg,'toi')
           tmpcfg=[];
           tmpcfg.toilim=cfg.toi;
           alldata = ft_redefinetrial(tmpcfg,alldata);
        else
           % nothing
        end
        %% 4.1.a Freqanalysis to get CSD
% TODO: ADD switch to mtmconvol 
        switch csdmethod
            case{'fft'}
                tmpcfg=[];
                tmpcfg.method='mtmfft';
                tmpcfg.output='powandcsd';
                tmpcfg.taper=taper;
                tmpcfg.keeptrials='no';
                tmpcfg.foi=foi;
                tmpcfg.tapsmofrq=tapsmofrq;

                csd_all = ft_freqanalysis(tmpcfg,alldata); % All data
                
            case{'tfr'}
                tmpcfg=[];
                tmpcfg.method='mtmconvol';
                tmpcfg.output='powandcsd';
                tmpcfg.taper=taper;
                tmpcfg.keeptrials='no';
                tmpcfg.foi=foi;
                tmpcfg.tapsmofrq=tapsmofrq;
                tmpcfg.toi=latency;
                tmpcfg.t_ftimwin=t_ftimwin;

                csd_all = ft_freqanalysis(tmpcfg,alldata); % All data
                csd_all.dimord='chan_freq';
                csd_all = rmfield(csd_all,'time');
        end
        %% 4.2.a Build Spatial Filter

        tmpcfg=[];
        tmpcfg.method='dics';
        tmpcfg.frequency=foi;
        tmpcfg.grid=lf2;
        tmpcfg.vol=vol;
        tmpcfg.dics.keepfilter='yes';
        tmpcfg.dics.realfilter='yes';
        tmpcfg.dics.lambda=lambda;
        tmpcfg.latency=latency;
        tmpcfg.dics.fixedori='yes';
        
        if ~isfield(csd_all,'grad') % If this is NOT a MEG structure
            tmpcfg.elec = elec; % Set the stuff you specified as "elec"
        end
            
        source = ft_sourceanalysis(tmpcfg,csd_all);

    case{'lcmv'}
        disp('Building the LCMV Filter')
        %% 4.1.b Timelockanalysis to get Covariance
        tmpcfg=[];
        tmpcfg.covariance = 'yes';
        tmpcfg.covariancewindow = covwin;
        
        avg = ft_timelockanalysis(tmpcfg,alldata);
        
        %% 4.2.b Build Spatial Filter
        
        tmpcfg=[];
        tmpcfg.method='lcmv';
        tmpcfg.grid=lf2;
        tmpcfg.vol=vol;
        tmpcfg.lcmv.keepfilter='yes';
        tmpcfg.lcmv.realfilter='yes';
        tmpcfg.lcmv.lambda=lambda;
        tmpcfg.lcmv.fixedori='yes';

        if ~isfield(avg,'grad') % If this is NOT a MEG structure
            tmpcfg.elec = elec; % Set the stuff you specified as "elec"
        end
        
        source = ft_sourceanalysis(tmpcfg,avg);
end


%% 5. Apply Spatial Filters to Data


virt=cell(1,length(indata.trial)); % Virtual Electrode

    for p=1:length(lf2.inside)    %loop grid points
    fprintf('Projecting single trials of Voxel %i\n', p)
        for k = 1:length(indata.trial)
            virt{k}(p,:)=source.avg.filter{lf2.inside(p)}*indata.trial{k};
        end %k
    end%i
out=virt;

  
%% 5.1. Make the fake FFT-Structure

    sp=[];
    
    % Make Labels
    if ~isempty(location)
        sp.label= cellstr(num2str(location));
    else
        for j = 1:length(lf2.inside);
            sp.label{j} = sprintf('%d',lf2.inside(j));
        end
    end

    % Plug in Data
    for t = 1:length(indata.trial)
            sp.time{t}=indata.time{1}; % Get the original Time back
    end

    sp.sourcedummy = source;
    sp.label=sp.label';
    sp.fsample=indata.fsample;
    sp.trial=out; %% trial-structre unneccessary
    for t=1:length(sp.trial)
        sp.time{t}=intime; % restore old time vector
    end

    % Add "electrode" strcture
    sp.elec.unit = 'mm';
    sp.elec.type = 'eeg';
    sp.elec.chanpos = source.pos(source.inside,:);
    sp.elec.elecpos = sp.elec.chanpos;
    sp.elec.label = sp.label; 
end
