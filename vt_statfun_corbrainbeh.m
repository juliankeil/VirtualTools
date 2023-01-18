function [s,cfg]=statfun_corbrainbeh(cfg, dat, design)
%This functions calculates a voxel-/sensorwise correlation of power with a
%vector of values (e.g. RT, questionnaire data etc.) passed on via the design matrix.
%(row of ivar). The length of the vector should correspond to the amount of data
%structures passed on by the statistic.
%
%Attention: the correlations are not calculated if there are NaNs in the
%vector containg the behavioural data. Remove according subjects for
%correct calculations.
%
% The external interface of this function has to be
%   [s,cfg] = statfun_corbrainbeh(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (iv) and the unit-of-observation (UO)
%          factor,  Nfac x Nreplications
%          IMPORTANT: the row with the ivar should contain the behavioural
%          data (e.g. RT) corresponding to each uvar. For biserial correlation,
%	   the ivar should only contain ones and zeros.
%
%The OUTPUT s contains the contains following fields:
%   -s.cor: Pearson Product Moment correlations / Point biserial correlation coefficients
%   -s.stat: t-values (formula taken from cor.test.R)
%   -s.prob: p-values (formula taken from cor.test.R)
%
% Configuration options:
%   cfg.dichotomous    = 'yes' or 'no', use point-biserial correlation for dichotomous variables (default='no')
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail = -1, 0, or 1, left, two-sided, or right (default=0)
%              cfg.tail in combination with cfg.computecritval='yes'
%              determines whether the critical value is computed at
%              quantile cfg.alpha (with cfg.tail=-1), at quantiles
%              cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%              quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification (see above):
%   cfg.ivar        = row number of the design that contains the labels of the conditions that must be
%                        compared (default=1). The labels are the numbers 1 and 2.
%   cfg.uvar        = row number of design that contains the labels of the UOs (subjects or trials)
%                        (default=2). The labels are assumed to be integers ranging from 1 to
%                        the number of UOs.
%
% Nathan Weisz, August 2008
% Mathis Kaiser, February 2015: added option for point-biserial correlation - code inspired by 
% http://www.mathworks.com/matlabcentral/fileexchange/11222-point-biserial-correlation/content/pointbiserial.m

%%
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=0;                end;
if ~isfield(cfg, 'dichotomous'),       cfg.dichotomous='no';      end;

%%

switch cfg.dichotomous
    case 'yes'
        if strcmp(cfg.computestat,'yes')
            behvec=design(cfg.ivar,:);
            if length(behvec) ~= size(dat,2)
                error('Mismatch between behvec and units of observation.');
            end
            sel0 = find(design(cfg.ivar,:)==0);
            sel1 = find(design(cfg.ivar,:)==1);
            n0  = length(sel0);
            n1  = length(sel1);
            if any(isnan(dat(:)))
                m0 = nanmean(dat(:,sel0), 2);
                m1 = nanmean(dat(:,sel1),2);
                stanDev = nanstd(dat, 1, 2); % calculated by dividing by n
            else
                m0 = mean(dat(:,sel0), 2);
                m1 = mean(dat(:,sel1),2);
                stanDev = std(dat, 1, 2);
            end
            s.cor = ((m1-m0)./stanDev)*sqrt(n0*n1/length(design)^2);
            for i=1:length(s.cor)
                if ~isnan(s.cor(i))
                    s.stat(i) = s.cor(i)/(sqrt((1-s.cor(i)^2)/(size(dat,2)-2)));
                else
                    s.stat(i) = NaN;
                end
            end
        end
%%
    case 'no'
        if strcmp(cfg.computestat,'yes')
            behvec=design(cfg.ivar,:);
            if length(behvec) ~= size(dat,2)
                error('Mismatch between behvec and units of observation.');
            end
            
            % Implement own correlation to avoid for-loops
            % TODO: behvec without variance results in zero-correlations,
            % should be NaNs
            dat=zscore(dat')';
            behvec=zscore(behvec);
            behmat=repmat(behvec,size(dat,1),1);
            
            df=size(dat,2)-1;
            
            s=[];
            s.cor=sum(dat.*behmat,2)/df;
            % Calculate t-values
            s.stat=(sqrt(repmat(df,length(s.cor),1)).*s.cor) ./ sqrt(1-s.cor.^2);
        end
end %switch
%%
if strcmp(cfg.computecritval,'yes')
    % also compute the critical values
    df      = size(dat,2)-2; %taken from cor.test
    s.df = df;
    if cfg.tail==-1
        s.critval = tinv(cfg.alpha,df);
    elseif  cfg.tail==0
        s.critval = [tinv(cfg.alpha/2,df),tinv(1-cfg.alpha/2,df)];
    elseif cfg.tail==1
        s.critval = tinv(1-cfg.alpha,df);
    end;
end

%%
if strcmp(cfg.computeprob,'yes')
    % also compute the p-values
    df      = size(dat,2)-2; %taken from cor.test
    s.df = df;
    if cfg.tail==-1
        s.prob = tcdf(s.stat,df);
    elseif  cfg.tail==0
        s.prob = 2*tcdf(-abs(s.stat),df);
    elseif cfg.tail==1
        s.prob = 1-tcdf(s.stat,df);
    end;
end

