%% Virtual Tools ANOVA Testskript
%% 1. Simulate some data
erp1_1 = [];
erp1_2 = [];
erp2_1 = [];
erp2_2 = [];
erp3_1 = [];
erp3_2 = [];
load Berlin_EEG_Head

%% Lets loop 20 Subjects
for vp = 1:20;
    %% 1. Create some data
    % 1.1. Make the ERP
    y=ones(1,20);
    g=gausswin(20);
    y2=y.*g';
    y3=y.*-g';
    begvec=zeros(1,198);
    endvec=zeros(1,98);
    y4=[begvec y2(1:end-1) y3(2:end)+.1 endvec]; % ERP with 2 Peaks and no Noise
    y5 = y4.*1.2; % Slightly larger ERP
    y6 = y4.*1.5; % Slightly larger ERP
    y7 = y4.*1.9; % Slightly larger ERP

    % 1.2. Make Pink Noise and add to ERP

    for t=1:50
        Nx=100000;
        B=[0.049922035 -0.95993537 0.050612699 -0.004408786]; % Whatever
        A=[1 -2.494956002 2.017265875 -0.522189400]; 
        nT60=round (log(1000)/1-max(abs(roots(A))));
        v=0+0.001.*randn(1,Nx+nT60);
        x=filter(B,A,v);
        x=x(nT60+1:end);
        x=resample(x,1,300);
        sig{t}=x;
    end

    %% 2. Create Simulated Recording

    %% 2.1. Scaling Factor for the ERP as it shouldn't look identical on all
    % channels, to see the order of the channels, look at the 1020.lay-file

    scalef=rand(1,length(elec_126.label));
    label=elec_126.label;
    %% 2.2. Make the data

    for t=1:50;
        for c=1:length(elec_126.label)
            data1_1.trial{t}(c,:)=sig{t}+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data1_2.trial{t}(c,:)=sig{t}+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data2_1.trial{t}(c,:)=sig{t}+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data2_2.trial{t}(c,:)=sig{t}+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data3_1.trial{t}(c,:)=sig{t}+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data3_2.trial{t}(c,:)=sig{t}+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
        end
        for c=[1,2,3,4,5,36,37,10,24,91,52];
            data1_1.trial{t}(c,:)=sig{t}+(y4*scalef(c))+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data1_2.trial{t}(c,:)=sig{t}+(y5*scalef(c))+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data2_1.trial{t}(c,:)=sig{t}+(y6*scalef(c))+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data2_2.trial{t}(c,:)=sig{t}+(y7*scalef(c))+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data3_1.trial{t}(c,:)=sig{t}+(((y6+y4)/2)*scalef(c))+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
            data3_2.trial{t}(c,:)=sig{t}+(((y7+y5)/2.1)*scalef(c))+(0+0.2.*randn(1,length(sig{1}))); % Signal is pink noise plus the scaled ERP plus some noise
        end
        data1_1.time{t}= -0.5:0.003:0.499;
        data1_2.time{t}= -0.5:0.003:0.499;
        data2_1.time{t}= -0.5:0.003:0.499;
        data2_2.time{t}= -0.5:0.003:0.499;
        data3_1.time{t}= -0.5:0.003:0.499;
        data3_2.time{t}= -0.5:0.003:0.499;
    end

    %% 2.4. Add the rest
    data1_1.label=elec_126.label;
    data1_1.fsample=length(sig{1});

    data1_2.label=elec_126.label;
    data1_2.fsample=length(sig{1});
    
    data2_1.label=elec_126.label;
    data2_1.fsample=length(sig{1});

    data2_2.label=elec_126.label;
    data2_2.fsample=length(sig{1});
    
    data3_1.label=elec_126.label;
    data3_1.fsample=length(sig{1});

    data3_2.label=elec_126.label;
    data3_2.fsample=length(sig{1});

    %% ERP
    erp1_1{vp}=ft_timelockanalysis([],data1_1);
    erp1_2{vp}=ft_timelockanalysis([],data1_2);
    
    erp2_1{vp}=ft_timelockanalysis([],data2_1);
    erp2_2{vp}=ft_timelockanalysis([],data2_2);
    
    erp3_1{vp}=ft_timelockanalysis([],data2_1);
    erp3_2{vp}=ft_timelockanalysis([],data2_2);
    
    %% TFR
    cfg=[];
    cfg.method='wavelet';
    cfg.output='pow';
    cfg.foi=[5:2:40];
    cfg.toi=[-.5:.01:.5];
    cfg.pad='nextpow2';

    tfr1_1{vp}=ft_freqanalysis(cfg,data1_1);
    tfr1_2{vp}=ft_freqanalysis(cfg,data1_2);
    
    tfr2_1{vp}=ft_freqanalysis(cfg,data2_1);
    tfr2_2{vp}=ft_freqanalysis(cfg,data2_2);
    
    tfr3_1{vp}=ft_freqanalysis(cfg,data2_1);
    tfr3_2{vp}=ft_freqanalysis(cfg,data2_2);
end

%% plot one sibject
plot(erp1_1{1}.time,erp1_1{1}.avg);

%% Prepare Neighbours for Stats

cfg=[];
cfg.elec = elec_126;
cfg.method = 'distance';
cfg.neighbourdist = 35;

neigh = ft_prepare_neighbours(cfg);

%% Compute Stats
cfg=[];
cfg.latency = [0.07 .16];
cfg.nIV1 = 3; % 3 Levels on the first factor
cfg.nIV2 = 2; % 2 Levels on the second factor
cfg.parameter = 'avg';
cfg.alpha = .05;
cfg.neighbours = neigh;
cfg.bf = 'yes';
cfg.correctm = 'cluster';
cfg.numrandomization = 3;
cfg.minnb = 2;

stats = vt_time_rmANOVA(cfg,erp1_1{:},erp1_2{:},erp2_1{:},erp2_1{:},erp3_1{:},erp3_1{:});

cfg.parameter = 'powspctrm';
cfg.frequency = [25 28];
stats = vt_freq_rmANOVA(cfg,tfr1_1{:},tfr1_2{:},tfr2_1{:},tfr2_1{:},tfr3_1{:},tfr3_1{:});
% Now we have F-Values, p-values (in the prob fields) and a neighbourhood-corrected mask (or FDR-corrected mask, if cfg.correctm = 'fdr')



%     %% HERE'S THE IMPORTANT PART!
% % Combining the Neighbour-Mask 
% for x = 1:length(stats.maskint(:))
%     if stats.maskint(x) == 0;
%         stats.maskint(x) = NaN;
%     end
% end

%%
figure;imagesc(stats.time,[],stats.statint)
figure;imagesc(stats.time,[],stats.statint.*stats.maskint)


stats.statint2 = stats.statint.*stats.maskint;
stats.probint2 = stats.probint.*stats.maskint;
%% Here's the Weathermap
% cfg.parameter are the values to be plotted
% cfg.probparameter will be used to compute the time-stability mask
% cfg.maskparameter will be multiplied with the time time-stability mask for additional correction.
cfg = [];
cfg.parameter = 'statint';
cfg.computemask = 'yes';
cfg.numconseq = 10;
cfg.alpha = 0.05;
cfg.probparameter = 'probint'; % pre-computed p-values 
cfg.maskparameter = 'maskint'; % pre-computed mask (e.g. neighborhood corrected)

vt_plot_weathermap(cfg,stats);

%%
GA1_1 = ft_timelockgrandaverage([],erp1_1{:});
GA1_2 = ft_timelockgrandaverage([],erp1_2{:});
GA2_1 = ft_timelockgrandaverage([],erp2_1{:});
GA2_2 = ft_timelockgrandaverage([],erp2_2{:});

cfg = [];
cfg.channel = 1:3;
cfg.layout = lay_126;

ft_singleplotER(cfg,GA1_1,GA1_2,GA2_1,GA2_2)
