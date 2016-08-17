%% Load Volume, MRI and Electrodes
%inpath = ('/home/keil/data/common/FieldTrip Tutorial/LabTests/Checker');
% load('Berlin_EEG_Head.mat');%, 'mri', 'elec_126');

cd /home/keil/data/projects/PTB/testpilot_6/

ptb = ft_read_mri('multi_pilot6-0004-0001.dcm');

%% Check

cfg=[];
cfg.method='ortho';

ft_sourceplot(cfg,ptb);

%% Realign

cfg=[];
cfg.method = 'interactive';
cfg.clim = [0 800];
ptb_r = ft_volumerealign(cfg,ptb);

%%

cfg=[];
cfg.output= {'brain' 'scalp' 'skull'};
mri_s = ft_volumesegment(cfg,ptb_r);

%%
cfg=[];
%cfg.tissue = 'brain';
cfg.interactive='no';
cfg.numvertices=400;

bnd_s = ft_prepare_mesh(cfg,mri_s); 

%%
figure;
triplot(bnd_s(3).pnt, bnd_s(3).tri,  [], 'edges');

%%

cfg=[];
cfg.method='dipoli';
%cfg.tissue='brain';
%cfg.hdmfile=vol;
cfg.conductivity = [0.3300 0.0041 0.3300];

vol_s = ft_prepare_headmodel(cfg,bnd_s);

% %% Make a slightly smaller Volume
% 
% vol2=vol_s;
% vol2.bnd(1).pnt = vol2.bnd(1).pnt.*.99;
% vol2.bnd(2).pnt = vol2.bnd(2).pnt.*.99;
% vol2.bnd(3).pnt = vol2.bnd(3).pnt.*.2;
%% Make a first Leadfield
% elec_126.unit='cm';
% vol.unit='cm';
% bnd_s(3).unit = 'cm';
 
cfg=[];
cfg.elec=elec_126;
cfg.vol=vol_s;

%posvec = 1:3:length(vol2.bnd(3).pnt);
% cfg.grid.pos = vol2.bnd(3).pnt(posvec,:);
cfg.grid.inside = 1:length(bnd_s(3).pnt);
cfg.grid.pos = bnd_s(3).pnt;

lfx=ft_prepare_leadfield(cfg);


%% Take another look

figure;
plot3(lfx.pos(lfx.inside,1),lfx.pos(lfx.inside,2),lfx.pos(lfx.inside,3),'*','LineWidth',5);
triplot(bnd_s(3).pnt, bnd_s(3).tri,  [], 'edges');

%% And now also mal a Layout

for l = 1:length(lfx.inside)
elec_lf.label{l}=num2str(lfx.inside(l));
end
elec_lf.elecpos=lfx.pos(lfx.inside,:);
elec_lf.chanpos=lfx.pos(lfx.inside,:);
elec_lf.type='eeg';

cfg=[];
cfg.projection='polar';
cfg.elec=elec_lf;
lay_lf = ft_prepare_layout(cfg);

cfg=[];
cfg.layout=lay_lf;
ft_layoutplot(cfg);


%%

% save berlin_elec_lf_vol_mri_hippstyle_hires.mat elec_126 lfx vol mri lay_lf
save berlin_elec_lf_vol_mri_hippstyle_vol.mat elec_126 lfx vol mri lay_lf