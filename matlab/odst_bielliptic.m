function [dv, tof, nburns, plan] = odst_bielliptic(r1, r2, rb, di)
%ODST_BIELLIPTIC Bi-elliptic transfer r1 -> rb -> r2, plane change at rb.
C = odst_const(); mu = C.MU;
a1 = (r1 + rb)/2;  a2 = (r2 + rb)/2;
vv = @(r, a) sqrt(mu*(2./r - 1./a));
dv1 = abs(vv(r1, a1) - sqrt(mu/r1));
dv2 = odst_combined(vv(rb, a1), vv(rb, a2), di);
dv3 = abs(vv(r2, a2) - sqrt(mu/r2));
t1 = pi*sqrt(a1^3/mu);  t2 = pi*sqrt(a2^3/mu);
dv = dv1 + dv2 + dv3;
tof = t1 + t2;
nburns = 3;
plan = [0,       vv(r1, a1),   0,  r1, dv1;
        t1,      vv(rb, a2),   di, rb, dv2;
        t1 + t2, sqrt(mu/r2),  0,  r2, dv3];
end
