function [aurtype, mlat, mlt, flux, energy] = empirical(doy, dFdt, hemi, premodel)

    function wt = season_weight(doy)
        wt = zeros(1, 4);
        if doy < 79
            wt(2) = 1 - (79 - doy)/90;
            wt(1) = 1 - wt(2);
        elseif doy < 171
            wt(3) = 1 - (171 - doy)/92;
            wt(2) = 1 - wt(3);
        elseif doy < 263
            wt(4) = 1 - (263 - doy)/92;
            wt(3) = 1 - wt(4);
        elseif doy < 354
            wt(1) = 1 - (354 - doy)/91;
            wt(4) = 1 - wt(1);
        else
            wt(2) = 1 - (79 - (doy-365))/90;
            wt(1) = 1 - wt(2);
        end
    end

    season = {'winter', 'spring', 'summer', 'fall'};
    aurtype = {'diff', 'mono', 'wave'};

    nszn = length(season);
    ntype = length(aurtype);
    nmlt = 96;
    nmlat = 160;

    mlt = linspace(0, 24, nmlt+1);
    mlt = mlt(1: nmlt);
    mlat = linspace(50, 90, nmlat/2+1);
    mlat = mlat(1: nmlat/2);
    if hemi == 'S'
        mlat = -mlat;
    end

    prob_coeff_1 = zeros(nszn, ntype, nmlt, nmlat);
    prob_coeff_2 = zeros(nszn, ntype, nmlt, nmlat);
    energy_coeff_1 = zeros(nszn, ntype, nmlt, nmlat);
    energy_coeff_2 = zeros(nszn, ntype, nmlt, nmlat);
    number_coeff_1 = zeros(nszn, ntype, nmlt, nmlat);
    number_coeff_2 = zeros(nszn, ntype, nmlt, nmlat);

    for iszn = 1: nszn
        for itype = 1: ntype
            filename = [premodel '/' season{iszn} '_prob_b_' aurtype{itype} '.txt'];
            opts = detectImportOptions(filename);
            opts.DataLines = [1 nmlt*nmlat];
            prob_coeffs = readmatrix(filename, opts);
            energy_coeffs = readmatrix([premodel '/' season{iszn} '_' aurtype{itype} '.txt'], opts);
            number_coeffs = readmatrix([premodel '/' season{iszn} '_' aurtype{itype} '_n.txt'], opts);
            for imlt = 1: nmlt
                prob_coeff_1(iszn, itype, imlt, :) = prob_coeffs((imlt-1)*nmlat+1: imlt*nmlat, 1);
                prob_coeff_2(iszn, itype, imlt, :) = prob_coeffs((imlt-1)*nmlat+1: imlt*nmlat, 2);
                energy_coeff_1(iszn, itype, imlt, :) = energy_coeffs((imlt-1)*nmlat+1: imlt*nmlat, 3);
                energy_coeff_2(iszn, itype, imlt, :) = energy_coeffs((imlt-1)*nmlat+1: imlt*nmlat, 4);
                number_coeff_1(iszn, itype, imlt, :) = number_coeffs((imlt-1)*nmlat+1: imlt*nmlat, 3);
                number_coeff_2(iszn, itype, imlt, :) = number_coeffs((imlt-1)*nmlat+1: imlt*nmlat, 4);
            end
        end
    end

    nt = length(dFdt);
    dFdt_rep = shiftdim(repmat(dFdt(:), 1, nszn, ntype, nmlt, nmlat), 1);

    prob_effect = repmat(prob_coeff_1, 1, 1, 1, 1, nt) + repmat(prob_coeff_2, 1, 1, 1, 1, nt) .* dFdt_rep;
    prob_effect(prob_effect < 0) = 0;
    prob_effect(prob_effect > 1) = 1;

    energy_all = repmat(energy_coeff_1, 1, 1, 1, 1, nt) + repmat(energy_coeff_2, 1, 1, 1, 1, nt) .* dFdt_rep;
    energy_all(energy_all < 0) = 0;
    energy_all(energy_all > 10) = 0.5;
    energy_all(energy_all > 5) = 5;

    number_all = repmat(number_coeff_1, 1, 1, 1, 1, nt) + repmat(number_coeff_2, 1, 1, 1, 1, nt) .* dFdt_rep;
    number_all(number_all < 0) = 0;
    number_all(number_all > 2e10) = 0;
    number_all(number_all > 2e9) = 1e9;

    if hemi == 'S'
        doy = 365 - doy;
    end
    wt = season_weight(doy);

    energy = zeros(ntype, nmlt, nmlat, nt);
    number = zeros(ntype, nmlt, nmlat, nt);
    for iszn = 1: nszn
        energy = energy + wt(iszn) * squeeze(energy_all(iszn, :, :, :, :));
        number = number + wt(iszn) * squeeze(number_all(iszn, :, :, :, :));
    end

    if hemi == 'N'
        energy = energy(:, :, nmlat/2+1: nmlat, :);
        number = number(:, :, nmlat/2+1: nmlat, :);
    end

    if hemi == 'S'
        energy = energy(:, :, 1: nmlat/2, :);
        number = number(:, :, 1: nmlat/2, :);
    end

    flux = energy;
    energy = flux ./ number / 1.6e-9;
    energy(~isfinite(energy)) = 0;
end