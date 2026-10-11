function dv = odst_combined(va, vb, ang)
%ODST_COMBINED Magnitude of the impulse between two velocity vectors of
%   magnitude va and vb separated by angle ang [rad] (law of cosines).
dv = sqrt(max(va.^2 + vb.^2 - 2.*va.*vb.*cos(ang), 0));
end
