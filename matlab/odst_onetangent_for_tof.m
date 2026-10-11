function k = odst_onetangent_for_tof(r1, r2, di, tofTarget, kmax)
%ODST_ONETANGENT_FOR_TOF Smallest apsis factor k whose one-tangent transfer
%   meets tofTarget. TOF decreases and dv increases with k, so the smallest
%   admissible k is the minimum-dv fast transfer. Returns [] if k > kmax needed.
if nargin < 5, kmax = 50; end
k0 = 1 + 1e-9;
f = @(k) tofOf(k) - tofTarget;
if f(k0) <= 0
    k = k0;
elseif f(kmax) > 0
    k = [];
else
    k = fzero(f, [k0, kmax], optimset('TolX', 1e-10));
end

    function t = tofOf(k)
        [~, t] = odst_onetangent(r1, r2, k, di);
    end
end
