function [flux_interp, energy_interp] = interp_satellite(ut, flux, energy, time, max_interval)
    [nmlat, nmlt, ~] = size(ut);
    nt = length(time);
    flux_interp = NaN(nmlat, nmlt, nt);
    energy_interp = NaN(nmlat, nmlt, nt);

    for imlat = 1: nmlat
        for imlt = 1: nmlt
            x = squeeze(ut(imlat, imlt, :));
            y1 = squeeze(flux(imlat, imlt, :));
            y2 = squeeze(energy(imlat, imlt, :));
            valid = ~(isnan(x) | isnan(y1) | isnan(y2));
            x = x(valid);
            y1 = y1(valid);
            y2 = y2(valid);
            if length(x) >= 2
                [~, ix] = sort(x);
                x = x(ix);

                if min(diff(x)) < max_interval
                    flux_interp(imlat, imlt, :) = interp1(x, y1(ix), time);
                    energy_interp(imlat, imlt, :) = interp1(x, y2(ix), time);
                end
            end
        end
    end
end