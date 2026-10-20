function loss = odst_gravity_loss(plan, isp, tw)
%ODST_GRAVITY_LOSS Finite-burn gravity loss of an impulsive burn plan [km/s].
%   Robbins (1966): dv*(w*tau)^2/24 per burn, with w the circular rate at the
%   burn radius and tau the burn time at constant thrust T = tw*m0*g0.
%   plan columns 4-5 hold each burn's radius [km] and dv [km/s]. Matches a
%   fixed-attitude burn within ~5 % for tw >= 0.2 and overestimates a
%   velocity-tracking burn about 2x (python/validate.py, section 4).
C = odst_const();
ve = isp*C.G0; m = 1; loss = 0;
for k = 1:size(plan, 1)
    r = plan(k, 4); dv = plan(k, 5);
    tau = isp/tw*m*(1 - exp(-dv/ve));
    loss = loss + dv*(sqrt(C.MU/r^3)*tau)^2/24;
    m = m*exp(-dv/ve);
end
end
