doy = 51;
hemi = 'N';
time = 6: 1/6: 11;

filename_sat = {'../f16_20140220.nc', '../f17_20140220.nc', '../f18_20140220.nc'};
coverage = 1/4;
max_interval = 1;

filename_grnd = '../themis20140220.nc';
flux_thres = 15;
energy_thres = 5;

aurtype = 1;
f = readmatrix('../omni2.lst', 'FileType', 'text');
dFdt = interp1((f(:, 2) - doy) * 24 + f(:, 3), f(:, 6).^(4/3) .* sqrt(f(:, 4).^2 + f(:, 5).^2).^(2/3) .* abs(sin(atan2(f(:, 4), f(:, 5))/2)).^(8/3), time);
premodel = '../premodel';
mlt_emp = 0: 0.1: 23.9;
mlat_emp = 50: 89;

ratio = [1 1/3 1/3 1/9];
downsampling = 'random';

alphaRadius = 0.2;
lb = 1;
a = -pi/18: pi/180: pi/18;
b = -pi/18: pi/180: pi/18;
r = pi/18: pi/180: 4*pi/18;
lev_flux = [0.5 1];
lev_energy = [1 5];
KN = 10;

min_data_points = 50;

mlt_sim = 0: 0.1: 23.9;
mlat_sim = 50: 89;
[r_sim, t_sim] = meshgrid(pi/2-deg2rad(mlat_sim), mlt_sim*pi/12);
x_sim = r_sim .* cos(t_sim);
y_sim = r_sim .* sin(t_sim);
loc_sim = [x_sim(:) y_sim(:)];

delta = pi/180;
centerweight = 4.01;
overlap = delta * 4.5;
weight = 1;
normalization = false;
rho = 1;
derivation = false;

filename_output = 'auroral_model.nc';

mlat_sat = ncread(filename_sat{1}, 'mlat');
mlt_sat = ncread(filename_sat{1}, 'mlt');
ut_sat = ncread(filename_sat{1}, 'ut_n');
flux_sat = ncread(filename_sat{1}, 'flux_n');
energy_sat = ncread(filename_sat{1}, 'energy_n');
for ifile = 2: length(filename_sat)
    ut_sat = cat(3, ut_sat, ncread(filename_sat{ifile}, 'ut_n'));
    flux_sat = cat(3, flux_sat, ncread(filename_sat{ifile}, 'flux_n'));
    energy_sat = cat(3, energy_sat, ncread(filename_sat{ifile}, 'energy_n'));
end

sat = grid_satellite(mlat_sat, mlt_sat, ut_sat, flux_sat, energy_sat, time, coverage);

[interp_flux, interp_energy] = interp_satellite(ut_sat, flux_sat, energy_sat, time, max_interval);

time_grnd = double(ncread(filename_grnd, 'time')) / 3600;
mlat_grnd = ncread(filename_grnd, 'mlat');
mlt_grnd = ncread(filename_grnd, 'mlt');
flux_grnd = ncread(filename_grnd, 'flux');
energy_grnd = ncread(filename_grnd, 'energy');

flux_grnd(flux_grnd < flux_thres) = NaN;
energy_grnd(energy_grnd < energy_thres) = NaN;

[grnd_mlat, grnd_mlt, grnd_flux, grnd_energy] = ground(time_grnd, mlat_grnd, mlt_grnd, flux_grnd, energy_grnd, time);

[~, emp_mlat, emp_mlt, emp_flux, emp_energy] = empirical(doy, dFdt, hemi, premodel);

emp = preprocess_empirical(emp_mlat, emp_mlt, emp_flux, emp_energy, aurtype, mlat_emp, mlt_emp, mlat_sat, mlt_sat, interp_flux, grnd_mlat, grnd_mlt, grnd_flux, alphaRadius);

r_sat = pi/2 - deg2rad(mlat_sat);
t_sat = mlt_sat * pi/12;
x_sat = r_sat .* cos(t_sat);
y_sat = r_sat .* sin(t_sat);

nt = length(time);
nmlt_sim = length(mlt_sim);
nmlat_sim = length(mlat_sim);

nlev_flux = length(lev_flux);
nlev_energy = length(lev_energy);
label_flux = linspace(0, 1, nlev_flux+1);
label_energy = linspace(0, 1, nlev_energy+1);
prob_flux = zeros(nt, nmlt_sim*nmlat_sim);
prob_energy = zeros(nt, nmlt_sim*nmlat_sim);
for it = 1: nt
    flux_interp = interp_flux(:, :, it);
    valid = flux_interp > 0;
    flux_interp = flux_interp(valid);
    label = zeros(1, length(flux_interp));
    for ilev = 1: nlev_flux-1
        label(flux_interp >= lev_flux(ilev) & flux_interp < lev_flux(ilev+1)) = label_flux(ilev+1);
    end
    label(flux_interp >= lev_flux(nlev_flux)) = 1;
    [~, score, ~] = predict(fitcknn([x_sat(valid) y_sat(valid)], label, 'NumNeighbors', KN), loc_sim);
    prob_flux(it, :) = score * label_flux';

    energy_interp = interp_energy(:, :, it);
    valid = energy_interp > 0;
    energy_interp = energy_interp(valid);
    label = zeros(1, length(energy_interp));
    for ilev = 1: nlev_flux-1
        label(energy_interp >= lev_energy(ilev) & energy_interp < lev_energy(ilev+1)) = label_energy(ilev+1);
    end
    label(energy_interp >= lev_energy(nlev_energy)) = 1;
    [~, score, ~] = predict(fitcknn([x_sat(valid) y_sat(valid)], label, 'NumNeighbors', KN), loc_sim);
    prob_energy(it, :) = score * label_energy';
end

[am, bm, rmin, rmax] = auroral_boundary(x_sat, y_sat, interp_flux, lb, a, b, r);

basis = setup_basis(am, bm, rmin, rmax, delta, centerweight, overlap, weight);

flux_sim = zeros(nt, nmlt_sim, nmlat_sim);
energy_sim = zeros(nt, nmlt_sim, nmlat_sim);

for it = 1: nt
    disp(it)
    mlt = {double(sat{it}.mlt), double(mlt_sat), grnd_mlt(:, :, it), emp{it}.mlt};
    mlat = {double(sat{it}.mlat), double(mlat_sat), grnd_mlat(:, :, it), emp{it}.mlat};
    flux = {double(sat{it}.flux), interp_flux(:, :, it), grnd_flux(:, :, it), emp{it}.flux};
    energy = {double(sat{it}.energy), interp_energy(:, :, it), grnd_energy(:, :, it), emp{it}.energy};

    n1 = length(mlt{1});
    for j = 2: 4
        nj = length(mlt{j});
        len = min(nj, round(n1 * ratio(j)));

        switch downsampling
        case 'sequential'
            idx = round(linspace(1, nj, len));
        case 'random'
            idx = randsample(nj, len);
        end
        mlt{j} = mlt{j}(idx);
        mlat{j} = mlat{j}(idx);
        flux{j} = flux{j}(idx);
        energy{j} = energy{j}(idx);
    end

    for j = 1: 4
        valid = flux{j} > 0 & energy{j} > 0;
        mlt{j} = mlt{j}(valid);
        mlat{j} = mlat{j}(valid);
        flux{j} = flux{j}(valid);
        energy{j} = energy{j}(valid);
    end

    loc = [cat(1, mlt{:}) cat(1, mlat{:})];
    flux = cat(1, flux{:});
    energy = cat(1, energy{:});

    [~, valid, ~] = unique(loc, 'rows');
    loc = loc(valid, :);
    flux = flux(valid);
    energy = energy(valid);

    if size(loc, 1) <= min_data_points
        continue
    end

    r = pi/2 - deg2rad(loc(:, 2));
    t = loc(:, 1) * pi/12;
    loc = [r.*cos(t) r.*sin(t)];
    flux = struct('loc', double(loc), 'val', double(log(flux)), 'err', ones(1, length(flux))');
    energy = struct('loc', double(loc), 'val', double(energy), 'err', ones(1, length(energy))');

    [y, W, Z, Q, phi] = constants(flux, basis{it}, normalization, rho, derivation);
    lambda = optimize(y, W, Z, Q, phi, exp(-9), exp(5), 5e-3);
    [d, c, rhoMLE, ~, M] = kriging(lambda, y, W, Z, Q, phi);
    [m, ~] = prediction(loc_sim, basis{it}, normalization, rho, lambda, Z, Q, phi, M, d, c, rhoMLE);
    flux_sim(it, :, :) = reshape(exp(m).*prob_flux(it,:)', nmlt_sim, nmlat_sim);

    [y, W, Z, Q, phi] = constants(energy, basis{it}, normalization, rho, derivation);
    lambda = optimize(y, W, Z, Q, phi, exp(-9), exp(5), 5e-3);
    [d, c, rhoMLE, ~, M] = kriging(lambda, y, W, Z, Q, phi);
    [m, ~] = prediction(loc_sim, basis{it}, normalization, rho, lambda, Z, Q, phi, M, d, c, rhoMLE);
    energy_sim(it, :, :) = reshape(m.*prob_energy(it,:)', nmlt_sim, nmlat_sim);
end

nccreate(filename_output, 'time', 'Dimensions', {'time', nt}, 'Datatype', 'single', 'Format', 'netcdf4')
nccreate(filename_output, 'mlt', 'Dimensions', {'mlt', nmlt_sim}, 'Datatype', 'single')
nccreate(filename_output, 'mlat', 'Dimensions', {'mlat', nmlat_sim}, 'Datatype', 'single')
nccreate(filename_output, 'flux', 'Dimensions', {'time', nt, 'mlt', nmlt_sim, 'mlat', nmlat_sim}, 'Datatype', 'single')
nccreate(filename_output, 'energy', 'Dimensions', {'time', nt, 'mlt', nmlt_sim, 'mlat', nmlat_sim}, 'Datatype', 'single')
ncwrite(filename_output, 'time', time)
ncwrite(filename_output, 'mlt', mlt_sim)
ncwrite(filename_output, 'mlat', mlat_sim)
ncwrite(filename_output, 'flux', flux_sim)
ncwrite(filename_output, 'energy', energy_sim)

function lambda = optimize(y, W, Z, Q, phi, xmin, xmax, tol)
    function likelihood = LK(l)
        [~, ~, ~, likelihood] = kriging(exp(l), y, W, Z, Q, phi);
        likelihood = -likelihood;
    end
    lambda = exp(fminbnd(@LK, log(xmin), log(xmax), optimset('FunValCheck', 'on', 'TolX', tol)));
end