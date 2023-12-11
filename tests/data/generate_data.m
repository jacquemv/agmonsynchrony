
function generate_data

T = icn_generate_poisson_si([0 200],[1 2],0.23,0.03);
save_timeseries('timeseries1.txt', T{1})
save_timeseries('timeseries2.txt', T{2})

T = icn_generate_poisson_si([0 500],[1 2],0.47,0.03);
save_timeseries('timeseries3.txt', T{1})
save_timeseries('timeseries4.txt', T{2})

%-----------------------------------------------------
function save_timeseries(fname, t)

fid = fopen (fname, "w");
fdisp (fid, t);
fclose (fid);
