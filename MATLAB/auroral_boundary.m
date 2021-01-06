function [am, bm, rmin, rmax] = auroral_boundary(x, y, flux, lb, a, b, r)
    [nx, ny] = size(x);

    na = length(a);
    nb = length(b);
    nr = length(r);

    r2 = r.^2;
    rlb2 = (r(1) - (r(2) - r(1))/2)^2;
    rub2 = (r(nr) + (r(nr) - r(nr-1))/2)^2;

    nt = size(flux, 3);
    am = zeros(1, nt);
    bm = zeros(1, nt);
    rmin = zeros(1, nt);
    rmax = zeros(1, nt);

    for it = 1: nt
        aurora = flux(:, :, it) >= lb;
        vote = zeros(na, nb, nr);

        for ix = 1: nx
            for iy = 1: ny
                if aurora(ix, iy)
                    for ia = 1: na
                        for ib = 1: nb
                            d2 = (x(ix, iy) - a(ia))^2 + (y(ix, iy) - b(ib))^2;
                            if d2 >= rlb2 && d2 <= rub2
                                [~, ir] = min(abs(r2 - d2));
                                vote(ia, ib, ir) = vote(ia, ib, ir) + 1;
                            end
                        end
                    end
                end
            end
        end

        [~, idx] = max(vote(:));
        [ia, ib, ~] = ind2sub([na nb nr], idx);
        am(it) = a(ia);
        bm(it) = b(ib);

        aur = flux(:, :, it) > 0;
        rm2 = (x(aur) - am(it)).^2 + (y(aur) - bm(it)).^2;
        rmin(it) = sqrt(min(rm2));
        rmax(it) = sqrt(max(rm2));
    end
end