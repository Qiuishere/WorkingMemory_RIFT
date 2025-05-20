function subjects = Start_up(rootdir)

addpath('/project/3018085.01/scripts/subfun')

addpath('/home/predatt/qiuhan/Documents/Toolboxs/fieldtrip-20250114')
addpath('/home/predatt/qiuhan/Documents/Toolboxs/fieldtrip-20250114/private')

addpath('/home/predatt/qiuhan/Documents/Toolboxs/MVPA-Light-master/startup')

subjects = datainfo(rootdir);
ft_defaults
ft_warning off
startup_MVPA_Light