addpath '/home/common/matlab/fieldtrip/qsub'
addpath '/project/3018085.01/scripts/PreProc'

for thesub = [6]
    req_mem   = 60*10^9;
    req_etime = 4*3600; % 2h
    jobs{thesub} = qsubfeval(@run_preproc_ica, thesub,  'memreq',  req_mem,  'timreq',  req_etime);
end