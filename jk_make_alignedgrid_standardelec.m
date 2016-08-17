%% Build a source grid based on the MNI brain, align it with the individual BEM and fit electrodes
% Based on the fieldtrip-tutorial
% Steps:
% 1. Load the MNI-Brain
% 1.1. Build a source grid on the MNI-Brain
% 2. Load the individual Brain
% 2.1. Bild an individual BEM-Model
% 2.2. Morph the template grid to the individual BEM
% 3. Load Polhemus Data
% 3.1. Fit Polhemus recorded locations to BEM
% 4. Build Leadfield

%% 1. Load the MNI brain
clear all
close all
clc
% Watch out! Adjust the folders
load('/home/keil/data/programs/Matlab//fieldtrip-20131028/template/headmodel/standard_mri.mat');

template = mri;
template.coordsys = 'neuromag'; % so that FieldTrip knows how to interpret the coordinate system
%% 1.1. Build a source grid
% segment the template brain and construct a volume conduction model (i.e. head model): this is needed
% for the inside/outside detection of voxels.
% cfg          = [];
% cfg.output={'brain' 'scalp' 'skull'};
% template_seg = ft_volumesegment(cfg, template);
% 
% % Create a Head-Model for the template (OR! Load the standard template)
% cfg          = [];
% cfg.method   = 'dipoli';
% template_vol = ft_prepare_headmodel(cfg, template_seg);
% template_vol = ft_convert_units(template_vol, 'cm'); % Convert the vol to cm, since the grid will also be expressed in cm

% If this fails, you can also load the standard template
load('/home/keil/data/programs/Matlab/fieldtrip-20131028/template/headmodel/standard_bem.mat');
template_vol = vol;
template_vol = ft_convert_units(template_vol, 'cm');
%% 1.2. Now we actually put the grid points to the vol
% construct the dipole grid in the template brain coordinates
% the source units are in cm
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -.2;
cfg.vol        = template_vol;
template_grid  = ft_prepare_sourcemodel(cfg);

%% 1.3. Plot to check
% make a figure with the template head model and dipole grid
figure
hold on
ft_plot_vol(template_vol, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% 2. Make individual Grid
% First we load our own MRI
% Again, check for folders

%cd '/home/keil/data/projects/PTB/MRI Complete/multisens_AF43784_19781002_multisens_AF43784_20131002/094615_t1_mpr_ns_sag_pat2_iso_asy_0004'

mri = ft_read_mri('JD45149-0004-0001.dcm'); % Paste name of first Scan here

%% 2.0.1.Brighten up the MRI
% The MRI from the PTB is super dark, we need to make this a bit brighter

for a = 1:size(mri.anatomy,1)
    for b = 1:size(mri.anatomy,2)
        for c = 1:size(mri.anatomy,3)
            if mri.anatomy(a,b,c) >=800
                mri.anatomy(a,b,c) = 800;
            end
        end
    end
end

%% 2.1. Build individual BEM-Model
% First Realign Volume. Set Nas LPA, RPA, and positive z-values

cfg=[];
cfg.method='interactive';
%cfg.clim=[0 800];
mri_r = ft_volumerealign(cfg,mri);
%% 2.1.1. Plot to check
cfg=[];
cfg.method='ortho';

ft_sourceplot(cfg,mri_r);

%% 2.1.2. Segment the individual volume
cfg=[];
cfg.output={'brain' 'scalp' 'skull'};
cfg.brainsmooth = 5;
cfg.scalpsmooth = 5;
cfg.brainthreshold = .5;
cfg.scalpthreshold = .05;
mri_s = ft_volumesegment(cfg,mri_r);
mri_s = ft_convert_units(mri_s, 'cm'); % Make sure its in cm!
%% 2.1.3. Plot to Check for holes

mri_p = mri_s;
mri_p.anatomy = mri_r.anatomy;
mri_p.brain = double(mri_p.brain);
mri_p.scalp = double(mri_p.scalp);
mri_p.skull = double(mri_p.skull);


cfg = [];
cfg.method='slice';
cfg.nslices = 30;
ft_sourceplot(cfg,mri_p); %only mri

cfg.funparameter = 'brain';
ft_sourceplot(cfg,mri_p); %segmented gray matter on top

cfg.funparameter = 'scalp';
ft_sourceplot(cfg,mri_p); %segmented white matter on top

cfg.funparameter = 'skull';
ft_sourceplot(cfg,mri_p); %segmented csf matter on top

%% 2.1.4. Create Headmodel

cfg=[];
cfg.interactive = 'no';
cfg.numvertices = 1000;

bnd_s = ft_prepare_mesh(cfg,mri_s);

cfg=[];
cfg.method='dipoli'; % Works on Mac!
%cfg.tissue='brain';
cfg.conductivity=[0.3300 0.0041 0.3300];

vol_s = ft_prepare_headmodel(cfg,bnd_s);

% 2.2. create the subject specific grid, using the template grid that has just been created
cfg = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri_r;
grid               = ft_prepare_sourcemodel(cfg);
 
%grid.pos = grid.pos*10; % Set right dimensions

%% All Done!
%
% Now:
% Use "grid" to build a lead-field and do source analysis
% Use "template_grid" to move the gridpoints back to the MNI-locations
% after the source analysis! Just copy the .pos xyz-grid and dim fields from the template to
% the source structure.
%% make a figure of the single subject headmodel, and grid positions
figure;
ft_plot_vol(vol_s, 'edgecolor', 'none'); alpha 0.4;
ft_plot_mesh(grid.pos(grid.inside,:));

%% 3. Load EEG-Positions:

load('Berlin_EEG_Head.mat','elec_128');
elec = elec_128;
elec.pnt = elec.chanpos;
elec = ft_convert_units(elec,'cm'); %WATCH OUT!
%% 3.0.1. Check the positions on the head surface. They don't look good, do they?

figure
triplot(vol_s.bnd(1).pnt, vol_s.bnd(1).tri,  [], 'edges');
hold
plot3(elec.pnt(:,1),elec.pnt(:,2),elec.pnt(:,3),'b*')
plot3(elec.pnt(1,1),elec.pnt(1,2),elec.pnt(1,3),'r*')
plot3(elec.pnt(2,1),elec.pnt(2,2),elec.pnt(2,3),'g*')
plot3(elec.pnt(3,1),elec.pnt(3,2),elec.pnt(3,3),'k*')

%% 3.1. Move the electrodes to match the head
% If you didn't make an error in the previous steps, the electrodes should
% already fit quite well.
% If they appear tilted by 90°, you set the Nasion, LPA and RPA instead of
% AC and PC
%
% If they appear tilted towards the back or front, you didn't quite hit AC
% and PC
%
% Ususally it should suffice to squish the sides from 1 to .95
% cfg=[];
% cfg.method='interactive';
% cfg.elec=elec;
% cfg.headshape=vol_s.bnd(1);
% 
% elec_new=ft_electroderealign(cfg);

cfg = [];
cfg.channel = 1:126;
cfg.plot = 'no';

elec_new = vt_eegrealign(cfg,vol_s,elec_128);
%% 3.1.1. Check the position
 figure
  triplot(vol_s.bnd(3).pnt, vol_s.bnd(3).tri,  [], 'faces_skin');
 hold
 triplot(vol_s.bnd(1).pnt, vol_s.bnd(1).tri,  [], 'edges');
 plot3(elec_new.elecpos(:,1),elec_new.elecpos(:,2),elec_new.elecpos(:,3),'r.','MarkerSize',20)
 camlight left
 view([-90 0]);
 
%% 4. Make the Leadfield.
elec_new.type='eeg';
 
cfg=[];
cfg.elec=elec_new;
cfg.vol=vol_s;
cfg.grid=grid;

lf=ft_prepare_leadfield(cfg);

%% 4.1. Plot to check

triplot(vol_s.bnd(1).pnt, vol_s.bnd(1).tri,  [], 'edges');
hold;
triplot(vol_s.bnd(3).pnt, vol_s.bnd(3).tri,  [], 'edges');
plot3(lf.pos(lf.inside,1),lf.pos(lf.inside,2),lf.pos(lf.inside,3),'*','LineWidth',5);
plot3(elec_new.elecpos(:,1),elec_new.elecpos(:,2),elec_new.elecpos(:,3),'r*');


%% And save!
% We'll need the tempalte grid later on, if we want to move all the grid
% points back to the original position for statistics and stuff.
cd '/home/keil/data/projects/PTB/Headmodels'
save([mri_r.hdr(1).Filename(3:9), '.mat'],'mri_r','mri_s','bnd_s','vol_s','grid','template_grid','elec_new','lf','-V7.3');