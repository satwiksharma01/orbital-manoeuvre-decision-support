function rb = odst_bielliptic_opt(r1, r2, di, rbmax, tofmax)
%ODST_BIELLIPTIC_OPT Optimal intermediate apoapsis rB for the bi-elliptic transfer.
%   min dv(rB)  s.t.  max(r1,r2) < rB <= rbmax,  TOF(rB) <= tofmax
%   dv(rB) is multi-modal, so a log-spaced grid search is refined with fminbnd.
%   TOF increases monotonically with rB, so the TOF constraint becomes an
%   upper bound on rB found with fzero. Returns [] if the TOF cannot be met.
if nargin < 5, tofmax = Inf; end
rlo = max(r1, r2)*(1 + 1e-9);
rhi = rbmax;
if tofOnly(rhi) > tofmax
    if tofOnly(rlo) > tofmax
        rb = [];
        return
    end
    rhi = fzero(@(x) tofOnly(x) - tofmax, [rlo, rhi]);
end
grid = logspace(log10(rlo), log10(rhi), 400);
dvs = arrayfun(@dvOnly, grid);
[dmin, i] = min(dvs);
lo = grid(max(i - 1, 1));
hi = grid(min(i + 1, numel(grid)));
rb = grid(i);
if hi > lo
    [x, f] = fminbnd(@dvOnly, lo, hi, optimset('TolX', 1e-3));
    if f <= dmin, rb = x; end
end

    function d = dvOnly(x)
        d = odst_bielliptic(r1, r2, x, di);
    end
    function t = tofOnly(x)
        [~, t] = odst_bielliptic(r1, r2, x, di);
    end
end
