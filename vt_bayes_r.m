function bf10 = vt_bayes_r(cfg);

% Computes Bayes Factors for Correlations
% Based on Wetzels & Wagenmaker 2012
% To interpret:
% bf10:
% > 100: Decisive Evidence for H1
% 30 - 100: Very Strong Evidence
% 10 - 30: Strong Evidence
% 3 - 10: Substantial Evidence
% 1 - 3: Anecdotal Evidece
% 1: No Evidence
% 0.33 - 1: Anecdodal Evidence for H0
% 0.1 - 0.33: Substantial Evidence
% 0.033 - 0.1: Strong Evidence
% 0.01 - 0.033: Very Strong Evidence
% < 0.01: Decisive Evidence for H0
%
% Use as: 
% cfg = [];
% cfg.r = Set Correlation Value
% cfg.n = Set Group Size
% cfg.g = Optional: Set Baysian Prior Value
% bf10 = vt_bayes_r(cfg);
%
% Julian Keil, 08.07.2015
%

%% Set Basics
r = cfg.r;
n = cfg.n;

%% Start Function
if isfield(cfg,'g')
    % If g is set, use that
    gval = g;
else % Else, compute g for a near infinite number and integrate
    for g = 1:100000
         int(g)=(1+g)^((n-2)/2)*(1+(1-r^2)*g)^(-(n-1)/2)*g^(-3/2)*exp(-n/(2*g));  
    end
    gval = trapz(int);
end

% Now Compute the Bayes Factor
bf10 = sqrt((n/2))/gamma(1/2)*gval;


end