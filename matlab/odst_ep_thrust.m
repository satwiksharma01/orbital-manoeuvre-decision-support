function F = odst_ep_thrust(m)
%ODST_EP_THRUST Electric thrust [N] from input power: F = 2 eta P / ve, as the jet
%   power is eta P (m.power_ep_kW, m.eff_ep, m.isp_ep). Mirrors ep_thrust() in
%   python/odst.py.
C = odst_const();
F = 2*m.eff_ep*m.power_ep_kW*1e3/(m.isp_ep*C.G0*1e3);
end
