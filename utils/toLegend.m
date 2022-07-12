% Set the legend for only some of the plotted items, specified by PlotNos,
% using the cell array Names, at location, Location. Returns a handle to
% the legend.
%
% Tom Okell, June 2022
% 
% h = toLegend(PlotNos,Names,Location)

function h = toLegend(PlotNos,Names,Location)

if nargin < 3; Location = 'NorthEast'; end

g = flip(get(gca,'children')); % Tweaked in Matlab2017a since sort no longer works!
h = legend(g(PlotNos),Names,'Location',Location);