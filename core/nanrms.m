function [rms, mu, sigma2] = nanrms(x,dim)
%NANRMS Root Mean Squared value, ignoring NaNs.
if nargin == 1 % let mean figure out which dimension to work along
mu = nanmean(x);
sigma2 = nanvar(x);
rms = sqrt(mu.*conj(mu) + sigma2);
else
mu = nanmean(x,dim);
sigma2 = nanvar(x,0,dim);
rms = sqrt(mu.*conj(mu) + sigma2);
end