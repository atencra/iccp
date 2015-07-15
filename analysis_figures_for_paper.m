analysis_figures_for_paper






The following are for data where the pair of neurons exhibited 
significant synchrony within 1 ms of 0 delay.

For fra parameter plots:

Move to directory D:\ICC_Pairs_Data\fra
then run:

fradata = iccpairs_ccpairs_select_max_plot_fra_params;

For population plots:

iccpairs_plot_fradata_cf_resptype(fradata);
iccpairs_plot_fradata_bw(fradata);
iccpairs_plot_fradata_q(fradata);



For STRF parameter plots:

[rfdata] = iccpairs_ccpairs_select_max_plot_strf_params(rfdata)
iccpairs_plot_rfdata_fr_pli(rfdata)
iccpairs_plot_ptdata_cf_q_latency(ptdata) 

[mtfdata] = iccpairs_ccpairs_select_max_strf_bmf(rfdata)
iccpairs_plot_mtfdata_btmf_bsmf(mtfdata)
iccpairs_plot_fiodata(fiodata)








