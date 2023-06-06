function h = prettyHxg(ax, data, displayname, color, edges, nbins)
    if nargin < 3
        displayname = 'data';
    end
    if nargin < 4
        color = 'k';
    end
    if (nargin < 6 && nargin < 5) || (isempty(edges) && isempty(nbins))
	    h = histogram(ax, data, 'displaystyle', 'stairs', 'normalization', 'probability', 'linewidth', 3, 'edgecolor', color, 'displayname', displayname);
    elseif nargin < 6 && ~isempty(edges)
		h = histogram(ax, data, 'displaystyle', 'stairs', 'normalization', 'probability', 'linewidth', 3, 'edgecolor', color, 'displayname', displayname, 'binedges', edges);
	elseif isempty(edges)
		h = histogram(ax, data, 'displaystyle', 'stairs', 'normalization', 'probability', 'linewidth', 3, 'edgecolor', color, 'displayname', displayname, 'numbins', nbins);
	end
end