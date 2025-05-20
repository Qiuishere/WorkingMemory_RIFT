
subj_id = 1;
addpath('/home/predatt/qiuhan/MEG/scripts/subfun')
subjects = Start_up;

%load(fullfile(subjects(subj_id).results, '03_decoding','decode_angle'));
for thestim = 1:2
    for theview = 1:2
        result = stat_decoding_on_freq{thestim, theview};
        mv_plot_result(result)
    end
end