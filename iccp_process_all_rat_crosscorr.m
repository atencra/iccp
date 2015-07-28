function iccpairs_process_all_rat_crosscorr

library('iccpairs');

cd('E:\Behavior - Temporal Processing - 2013-2014\2013-11-25-trained-rat-tm20-acute\single_multi_unit_strf');
cc_batch_calc_chan_pair_crosscorr;


cd('E:\Behavior - Temporal Processing - 2013-2014\2013-12-19-trained-rat-tm20-acute\single_multi_unit_strf');
cc_batch_calc_chan_pair_crosscorr;


cd('E:\Behavior - Temporal Processing - 2013-2014\2014-3-6-control-rat\singe_multi_unit_strf');
cc_batch_calc_chan_pair_crosscorr;

cd('E:\Behavior - Temporal Processing - 2013-2014\2014-3-25-trained-rat-tm20-acute\single_multi_unit_strf');
cc_batch_calc_chan_pair_crosscorr;

exit;



