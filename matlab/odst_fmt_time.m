function s = odst_fmt_time(t)
%ODST_FMT_TIME Format a duration in seconds as hours (< 3 days) or days.
if t < 72*3600
    s = sprintf('%.2f h', t/3600);
else
    s = sprintf('%.1f d', t/86400);
end
end
