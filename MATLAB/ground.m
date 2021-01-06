function [mlatout, mltout, fluxout, energyout] = ground(timein, mlat, mlt, flux, energy, timeout)
    ntout = length(timeout);
    [nmlt, nmlat, ntin] = size(mlat);
    mlatout = zeros(nmlt, nmlat, ntout);
    mltout = zeros(nmlt, nmlat, ntout);
    fluxout = zeros(nmlt, nmlat, ntout);
    energyout = zeros(nmlt, nmlat, ntout);

    s = 1;
    while s <= ntout && timeout(s) < timein(1)
        s = s + 1;
    end

    j = 1;
    for i = s: ntout
        while j <= ntin-1 && timein(j+1) <= timeout(i)
            j = j + 1;
        end
        if j == ntin
            break
        end
        mlatout(:, :, i) = mlat(:, :, j);
        mltout(:, :, i) = mlt(:, :, j);
        fluxout(:, :, i) = flux(:, :, j);
        energyout(:, :, i) = energy(:, :, j);
    end
end