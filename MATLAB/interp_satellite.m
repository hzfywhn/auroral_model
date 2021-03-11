function [flux_interp, energy_interp] = interp_satellite(ut, flux, energy, time)
    [nmlat, nmlt, ~] = size(ut);
    nt = length(time);
    flux_interp = NaN(nmlat, nmlt, nt);
    energy_interp = NaN(nmlat, nmlt, nt);

    for imlat = 1: nmlat
        for imlt = 1: nmlt
            t = squeeze(ut(imlat, imlt, :));

            f = squeeze(flux(imlat, imlt, :));
            valid = ~(isnan(t) | isnan(f));
            x = t(valid);
            f = f(valid);
            [x, uniq, ~] = unique(x);
            if length(x) >= 2
                flux_interp(imlat, imlt, :) = interp1(x, f(uniq), time);
            end

            e = squeeze(energy(imlat, imlt, :));
            valid = ~(isnan(t) | isnan(e));
            x = t(valid);
            e = e(valid);
            [x, uniq, ~] = unique(x);
            if length(x) >= 2
                energy_interp(imlat, imlt, :) = interp1(x, e(uniq), time);
            end
        end
    end
end