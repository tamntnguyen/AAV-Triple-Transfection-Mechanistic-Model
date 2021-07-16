function [mu,c] = calc_growth_rate(t)
%% Cell growth fitted to cell growth curve as a function of time
g   = [49427.9581111241, 0.0112812472592155];
c   = 1e6 + g(1)/g(2)*(1 - exp(-g(2).*t));
mu  = g(1)*exp(-g(2).*t)./c;
