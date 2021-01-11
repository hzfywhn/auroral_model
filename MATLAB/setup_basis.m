function basis = setup_basis(am, bm, rmin, rmax, sponge, delta, centerweight, overlap, weight)
    nt = length(rmax);
    basis = cell(1, nt);

    for it = 1: nt
        inner = rmin(it) - sponge;
        if inner < 0
            inner = 0;
        end
        outer = rmax(it) + sponge;
        xy = -outer: delta: outer;
        nxy = length(xy);
        [x, y] = meshgrid(xy, xy);

        radius2 = x.^2 + y.^2;
        valid = radius2 >= inner^2 & radius2 <= outer^2;

        x = zeros(1, nxy^2);
        y = zeros(1, nxy^2);
        index = zeros(nxy, nxy);
        idx = 1;
        for i = 1: nxy
            for j = 1: nxy
                if valid(i, j)
                    x(idx) = xy(i);
                    y(idx) = xy(j);
                    index(i, j) = idx;
                    idx = idx + 1;
                end
            end
        end

        m = idx - 1;
        x = x(1: m);
        y = y(1: m);
        con = NaN(m, 4);

        idx = 1;
        for i = 1: nxy
            for j = 1: nxy
                if valid(i, j)
                    if i >= 2 && valid(i-1, j)
                        con(idx, 1) = index(i-1, j);
                    end
                    if j >= 2 && valid(i, j-1)
                        con(idx, 2) = index(i, j-1);
                    end
                    if i <= nxy-1 && valid(i+1, j)
                        con(idx, 3) = index(i+1, j);
                    end
                    if j <= nxy-1 && valid(i, j+1)
                        con(idx, 4) = index(i, j+1);
                    end
                    idx = idx + 1;
                end
            end
        end

        basis{it} = {struct('loc', [x'+am(it) y'+bm(it)], 'connect', con, 'centerweight', centerweight, 'delta', overlap, 'alpha', weight)};
    end
end