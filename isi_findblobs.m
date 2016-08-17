sourcedummy.avg.pow2 = sourcedummy.avg.pow;

sourcedummy.avg.pow2(sourcedummy.inside) = stats.stat(:,mfreq);

sourcedummy.avg.pow2(sourcedummy.outside) = 0;

%%

mat4clust = reshape(sourcedummy.avg.pow2, sourcedummy.dim);

%%

test=bwlabeln(mat4clust,26);
unique(test(:))

%%
roi=test==1;
length(find(roi(:)==1))

%%
index = [];

for r = 1:size(test,1)
    for c = 1:size(test,2)
        for d = 1:size(test,3)
            index(end+1) = sub2ind(size(test),r,c,d);
        end
    end
end

%%
clustermat = test(index)';

%%
sourcedummy.avg.pow3 = sourcedummy.avg.pow;
sourcedummy.avg.pow3(sourcedummy.inside) = stats.stat(:,mfreq);
sourcedummy.avg.pow3(sourcedummy.outside) = NaN;
sourcedummy.avg.pow3 = sourcedummy.avg.pow3.*(clustermat == 2);

%%

tmpcfg = [];
tmpcfg.parameter = 'avg.pow3';
tmpcfg.downsample = 1;

sourcedummy_i = ft_sourceinterpolate(tmpcfg,sourcedummy,mri);
    
    if isfield(sourcedummy_i,'time')
        sourcedummy_i = rmfield(sourcedummy_i,'time');
    end
    
%% And Plot the Surface
tmpcfg = [];
tmpcfg.method='ortho';
tmpcfg.surffile = 'surface_wminf_both.mat';
tmpcfg.funparameter = 'avg.pow3';
tmpcfg.funcolorlim = [0 20];
tmpcfg.maskparameter = tmpcfg.funparameter;
%tmpcfg.atlas = '/home/keil/data/programs/Matlab/fieldtrip-20131028/template/atlas/afni/TTatlas+tlrc.BRIK';%aal/ROI_MNI_V4.nii';
%tmpcfg.atlas = '/home/keil/data/programs/Matlab/fieldtrip-20131028/template/atlas/aal/ROI_MNI_V4.nii';
tmpcfg.coordsys = 'mni';
tmpcfg.camlight = 'off';

ft_sourceplot(tmpcfg,sourcedummy_i);

%%
barmask = NaN(990,1);
barmask(clustermat==2) = sourcedummy.avg.pow3(clustermat==2);
barmask = barmask(sourcedummy.inside);

[chans] = find(barmask > 0);
chans = unique(chans);

mfreq=10; % in case of beta or gamma
%[peakval peakchan] = max(mval);

figure;
cfg=[];
cfg.xlabel = 'visfix > visvar > tacfix > tacvar';
cfg.ylabel = sprintf('Power at %d Hz', round(stats.freq(mfreq)));
cfg.ylim = [-.2 .2];

jk_bar_err(cfg,mean(reshape(squeeze(stats.means(chans,mfreq,:,:,:)),length(chans),4)),...
    mean(reshape(squeeze(stats.sems(chans,mfreq,:,:,:)),length(chans),4)));
