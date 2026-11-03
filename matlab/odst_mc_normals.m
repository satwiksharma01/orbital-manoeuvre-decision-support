function Z = odst_mc_normals()
%ODST_MC_NORMALS The Monte Carlo runs: 4000 x 10 standard normals, Box-Muller on
%   consecutive uniform pairs from odst_mrg32k3a, run by run. Columns: the Isp error, the
%   thrust error, then the two pointing components of the first, second and third burn;
%   the fourth burn's (the nodal drift has four) come from the stream after the first
%   8 x 4000 uniforms. The seed is fixed, so every option, mission and engine uses the same
%   runs (mc_normals in odst.py).
persistent Zc
if isempty(Zc)
    N = 4000;
    u = odst_mrg32k3a(10*N);
    r = sqrt(-2*log(u(1:2:8*N))); a = 2*pi*u(2:2:8*N);
    r4 = sqrt(-2*log(u(8*N + 1:2:end))); a4 = 2*pi*u(8*N + 2:2:end);
    Zc = [reshape([r.*cos(a), r.*sin(a)].', 8, N).', r4.*cos(a4), r4.*sin(a4)];   % row k: the pairs 4k-3 .. 4k, cos then sin
end
Z = Zc;
end
