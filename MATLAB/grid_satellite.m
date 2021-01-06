function gridded = grid_satellite(mlat, mlt, ut, flux, energy, time, coverage)
    reps = size(ut, 3);
    mlat = repmat(mlat, 1, 1, reps);
    mlt = repmat(mlt, 1, 1, reps);

    nt = length(time);
    gridded = cell(1, nt);

    for i = 1: nt
        mask = ut > time(i) - coverage & ut < time(i) + coverage;
        gridded{i} = struct('mlat', mlat(mask), 'mlt', mlt(mask), 'flux', flux(mask), 'energy', energy(mask));
    end
end