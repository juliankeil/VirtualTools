function [pdv] = jk_phase_deviance(cfg,indata);

% computes single trial deviance from dircular difference in phase values
% i.e. takes the phase of two channels, computes the circular difference in
% phase angles and compares this to the mean phase difference over trials
% 
% Conceptually related to ITC
%
% input = fourier values
% output = difference values in rad (-pi to pi)
%
% cfg.freq = frequency
% cfg.trials = trials
% cfg.chancmb = Channel Combination for Deviance
%
%
% Ver 1, 02.09.2013
% Ver 2, 16.10.2013: Added Channel Combination option
% Ver 2.1. 23.10.2013: Fixed Bug for multiple Channels
% (c) Julian Keil, 2013

%% Set CFGs

if isfield(cfg,'freq')
    for f = 1:length(cfg.freq)
        freq(f) = nearest(indata.freq,cfg.freq(f));
    end
else
    freq = 1:length(indata.freq);
end

if isfield(cfg,'trials')
    trials = cfg.trials;
else
    trials = size(indata.fourierspctrm,1);
end

%% 0. Prepare Stuff

if isfield(cfg,'chancmb')
    for l = 1:size(cfg.chancmb,1)
        pairs = cfg.chancmb;
        labels{l} = [num2str(cfg.chancmb{l,1}),'_',num2str(cfg.chancmb{l,2})];
    end
    %labels = unique(labels);
    
else
    chan = 1:length(indata.label);
    [idx2,idx1]=find(true(numel(chan),numel(chan)));
    pairs = [reshape(chan(idx1),[],1),reshape(chan(idx2),[],1)];

    tmplabels = [reshape(indata.label(idx1),[],1),reshape(indata.label(idx2),[],1)];

    for l = 1:length(tmplabels)
        lab1 = str2mat(tmplabels{l,1});
        lab2 = str2mat(tmplabels{l,2});
        labels{l} = cellstr([lab1,'_',lab2]);
    end
end

%% 1. Compute the Phase for all channels and trials

tmpphase = angle(indata.fourierspctrm);

%% 2. For each trial, compute the phase difference between all channels
clear tmpdiff

%tic
for f=1:length(freq)
    fprintf('Working on Frequency %4.2f \n',indata.freq(freq(f)))
    for t=1:trials
        fprintf('Working on Trial %d \n',t)
        for p=1:length(labels)
            %fprintf('Working on Pair %d \n',p)
            tmp = strcmpi(pairs(p,1),indata.label);
            p1 = find(tmp==1);
            tmp = strcmpi(pairs(p,2),indata.label);
            p2 = find(tmp==1);
            tmpdiff(t,p,f) = (circ_dist(tmpphase(t,p1,freq(f)),tmpphase(t,p2,freq(f))));
        end
    end
end
%toc

%% 3. For each combination, compute the mean over trials
clear tmpmean 

tmpmean = (circ_mean(tmpdiff));

%% 4. For each combination, compute the deviance from the mean
clear pdv outdist

for f=1:length(freq) % Frequency
    fprintf('Computing the Deviance on Frequency %4.2f \n',indata.freq(f))
    for t=1:size(tmpdiff,1) % Trials
        fprintf('Computing the Deviance on Trial %d \n',t)
        for p=1:length(labels) % Deviance for each trial from mean
            %fprintf('Working on Pair %d \n',p)
            tmp = strcmpi(pairs(p,1),indata.label);
            p1 = find(tmp==1);
            tmp = strcmpi(pairs(p,2),indata.label);
            p2 = find(tmp==1);
            outdist(t,p,f) = (circ_dist(tmpdiff(t,p,f),tmpmean(1,p,f)));
        end
    end
end

%% A little wrap-up

pdv.deviance = round(outdist.*1e+16)./1e+16; % Avoid nasty rounding error at 0
pdv.label = labels;
pdv.dimord = indata.dimord;
pdv.freq = indata.freq(freq);