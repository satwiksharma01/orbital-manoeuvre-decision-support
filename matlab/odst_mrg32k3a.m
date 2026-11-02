function u = odst_mrg32k3a(n, seed)
%ODST_MRG32K3A n uniforms in (0, 1) from L'Ecuyer's MRG32k3a (Operations Research 47, 1999).
%   In floating point like his reference code: every product is an exact integer below
%   2^53, so MATLAB, Python and the web draw identical numbers (mrg32k3a in odst.py).
%   seed: six integers, default 12345 each.
if nargin < 2, seed = 12345*ones(1, 6); end
m1 = 4294967087; m2 = 4294944443;
s1 = seed(1:3); s2 = seed(4:6);
u = zeros(n, 1);
for k = 1:n
    p1 = 1403580*s1(2) - 810728*s1(1); p1 = p1 - fix(p1/m1)*m1; if p1 < 0, p1 = p1 + m1; end
    p2 = 527612*s2(3) - 1370589*s2(1); p2 = p2 - fix(p2/m2)*m2; if p2 < 0, p2 = p2 + m2; end
    s1 = [s1(2) s1(3) p1]; s2 = [s2(2) s2(3) p2];
    if p1 > p2, u(k) = (p1 - p2)*2.328306549295727688e-10; else, u(k) = (p1 - p2 + m1)*2.328306549295727688e-10; end
end
end
