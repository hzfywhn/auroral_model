function emp_new = preprocess_empirical(emp_mlat, emp_mlt, emp_flux, emp_energy, aurtype, mlat_emp, mlt_emp, interp_mlat, interp_mlt, interp_flux, grnd_mlat, grnd_mlt, grnd_flux, alphaRadius)
    emp_flux = squeeze(emp_flux(aurtype, :, :, :));
    emp_energy = squeeze(emp_energy(aurtype, :, :, :));

    nmlt0 = length(emp_mlt);
    emp_mlt(nmlt0+1) = 24;
    emp_flux(nmlt0+1, :, :) = emp_flux(1, :, :);
    emp_energy(nmlt0+1, :, :) = emp_energy(1, :, :);

    nt = size(emp_flux, 3);
    nmlt = length(mlt_emp);
    nmlat = length(mlat_emp);
    [mlat_emp, mlt_emp] = meshgrid(mlat_emp, mlt_emp);

    flux = zeros(nmlt, nmlat, nt);
    energy = zeros(nmlt, nmlat, nt);
    for it = 1: nt
        flux(:, :, it) = interp2(emp_mlt, emp_mlat, emp_flux(:, :, it)', mlt_emp, mlat_emp);
        energy(:, :, it) = interp2(emp_mlt, emp_mlat, emp_energy(:, :, it)', mlt_emp, mlat_emp);
    end
    emp_flux = flux;
    emp_energy = energy;

    emp_mlt = mlt_emp;
    emp_mlat = mlat_emp;

    r = pi/2 - deg2rad(emp_mlat);
    t = emp_mlt * pi/12;
    x = r .* cos(t);
    y = r .* sin(t);

    emp_new = cell(1, nt);
    for it = 1: nt
        valid_interp = interp_flux(:, :, it) > 0;
        valid_grnd = grnd_flux(:, :, it) > 0;
        grnd_mlt_i = grnd_mlt(:, :, it);
        grnd_mlat_i = grnd_mlat(:, :, it);
        mlt_obs = [interp_mlt(valid_interp); grnd_mlt_i(valid_grnd)];
        mlat_obs = [interp_mlat(valid_interp); grnd_mlat_i(valid_grnd)];
        loc_obs = unique([mlt_obs mlat_obs], 'row');
        r = pi/2 - deg2rad(loc_obs(:, 2));
        t = loc_obs(:, 1) * pi/12;
        valid = ~inShape(alphaShape(double(r.*cos(t)), double(r.*sin(t)), alphaRadius), x, y);
        emp_flux_i = emp_flux(:, :, it);
        emp_energy_i = emp_energy(:, :, it);
        emp_new{it} = struct('mlt', atan2(y(valid), x(valid)) * 12/pi, 'mlat', 90 - rad2deg(sqrt(x(valid).^2 + y(valid).^2)), 'flux', emp_flux_i(valid), 'energy', emp_energy_i(valid));
    end
end