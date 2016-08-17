function [out] = vt_effect_size(cfg,F,df1,df2);

% This function calculates different effect sizes from ANOVAs
% The input requires a list of F-Values and degrees of freedom.
% Watch out! Lists of F and df need to be equally long. 
% Output can be - depending on cfg.method - partial eta square, partial
% omega squared or partial epsilon squared.
% All are scaled between 0 and 1, i.e.
% 0.1 = small effect size
% 0.3 = medium effect size
% 0.5 = large effect size
% Check en.wikipedia.org/wiki/Effect_size for more info.
%
% cfg.method = 'eta'; % the widely used but perhaps biased partial eta squared.
% cfg.method = 'omega'; % unbiased omega squared
% cfg.method = 'epsilon'; % unbiased epsilon squared
% cfg.method = 'all'; % compute all of the above
%
% Julian Keil, 2016
% Ver 1.: 20.04.2016: First implementation

%% Set the cfgs
if exist('cfg');
    if isfield(cfg,'method');
        if strcmpi(cfg.method,'eta')
            meth_flag = 1;
        elseif strcmpi(cfg.method,'omega')
            meth_flag = 2;
        elseif strcmpi(cfg.method,'epsilon')
            meth_flag = 3;
        elseif strcmpi(cfg.method,'all')
            meth_flag = 4;
        else
            meth_flag = 4;
        end % method
    else
        meth_flag = 4;
    end % field
else
    meth_flag = 4;
end % exist

%% Compute the Metric

if meth_flag == 1;
    for i = 1:length(F);
        out.eta(i) = (F(i)*df1(i))/((F(i)*df1(i))+df2(i));
    end
elseif meth_flag == 2;
    for i = 1:length(F);
        out.omega(i) = (F(i)-1)/(F(i)+((df2(i)+1)/df1(i)));
    end
elseif meth_flag == 3;
    for i = 1:length(F);
        out.epsilon(i) = (F(i)-1)/(F(i)+(df2(i)/df1(i)));
    end
elseif meth_flag ==4;
    for i = 1:length(F);
        out.eta(i) = (F(i)*df1(i))/((F(i)*df1(i))+df2(i));
        out.omega(i) = (F(i)-1)/(F(i)+((df2(i)+1)/df1(i)));
        out.epsilon(i) = (F(i)-1)/(F(i)+(df2(i)/df1(i)));
    end
end
