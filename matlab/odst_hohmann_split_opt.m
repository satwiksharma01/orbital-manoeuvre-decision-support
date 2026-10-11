function s = odst_hohmann_split_opt(r1, r2, di)
%ODST_HOHMANN_SPLIT_OPT Numerical optimisation of the plane-change split.
%   s* = argmin_{s in [0,1]} dv_Hohmann(s), s = fraction of di at the r1 burn.
if di <= 0 || abs(r2 - r1) < 1e-6
    s = double(r2 < r1);
    return
end
s = fminbnd(@(x) odst_hohmann(r1, r2, di, x), 0, 1, optimset('TolX', 1e-8));
end
