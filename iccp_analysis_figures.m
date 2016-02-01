%% ICC paired recordings analysis: Functions to make figures for paper
%% ICC Pairs Analysis
% The steps below were used to analyze central of the inferior colliculus (ICC) 
% paired neuron data. We recorded from the ICC using multi-channel electrodes, and from
% each electrode channel we extracted single units. In some cases, we found
% 2-3 single units from one channel. We then analyzed these pairs of
% neurons to determine how local circuits are functionally organized in the
% ICC.
%
% The ICC is important because it is the main processing center in the
% auditory midbrain. It is an obligatory processing center for all
% information in the ascending auditory system.
%
% For the experiments we presented pure tones, complex tones, broadband
% dynamic moving ripple stimuli. From the pure tones we estimated frequency
% response areas, and from the ripples we estimated spectrotemporal
% receptive fields.




%% Analysis Procedure
% # Raw data traces are spike-sorted using the SpikeSort1.3 program.
% # Save .spk and .mdl files in *-spk.mat files.
% # Get trigger for stimulus and save in *-trig.mat file.
% # Fix inter-spike intervals using: spknew = iccp_isi_filtered_spk(spk, delay); Or, if *-strfcmb-pairs.mat files already exist, then you can modify them using iccp_fix_strf_spk_isi.m
% # Estimate spectrotemporal receptive fields using calculate_strf.m
% # Estimate strf parameters using strf_parameters.m
% # Estimate modulation transfer function parameters for each strf using ripple_transfer_function.m
% # Estimate pure tone params using strf_pure_tone_params.m
% # Estimate nonlinearity using iccp_batch_sta_nonlinearity.m
% # Fit curves to nonlinearities using iccp_batch_sta_fio_fitParams.m
% # Estimate functional connectivity using iccp_batch_channel_corr.m 
% # This is a batch function which processes all data files in the current folder. 
% # It calls two functions: [pairstrains] = iccp_get_spk_paired_data(spk); [ccpairs] = iccp_calc_spk_crosscorr(pairstrains);
% # After this processing is complete, plotting function may be used to analyze the results.


%% Plotting functions:
% # [nb_nb_ne, ne_nb_ne] = iccp_plot_pairs_fra_resptype(nb_nb_ne, ne_nb_ne)



%% Synchrony analysis:
% # Get cross-corr data struct: iccp_batch_channel_corr.m
% # Add ccc to ccpairs: ccpairs = add_ccc_to_pairs_crosscorr(ccpairs)
% # Analyze pos/neg corr peaks: ccpairsPosNeg = iccp_ccpairs_cross_covariance_analysis(ccpairs)
% # Write code to combine all ccpairs struct into one long struct array: ccpairs = iccp_combine_ccpairs;
% # Plot cross-covariance parameters: iccp_plot_crosscorr_pd_hw_ccc_pos_neg(ccpairs, showstats)
% # The following are for data where the pair of neurons exhibited significant synchrony within 1 ms of 0 delay.
% # For fra parameter plots: Move to directory D:\ICC_Pairs_Data\fra then run: fradata = iccpairs_ccpairs_select_max_plot_fra_params;
% # For population plots:
% # iccpairs_plot_fradata_cf_resptype(fradata);
% # iccpairs_plot_fradata_bw(fradata);
% # iccpairs_plot_fradata_q(fradata);
% # data = iccp_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
% # data = iccp_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs)
% # data = iccp_plot_strf_bmf_crosscorr_pos_neg(ccpairs)
% # data = iccp_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs)

%% STRF parameter plots:
% * [rfdata] = iccpairs_ccpairs_select_max_plot_strf_params(rfdata)
% * iccpairs_plot_rfdata_fr_pli(rfdata)
% * iccpairs_plot_ptdata_cf_q_latency(ptdata) 
% * [mtfdata] = iccpairs_ccpairs_select_max_strf_bmf(rfdata)
% * iccpairs_plot_mtfdata_btmf_bsmf(mtfdata)
% * iccpairs_plot_fiodata(fiodata)

%% Plots with ccpairs struct array. Assumes that pos and neg temporal interactions have been analyzed. See Synchony Analysis above.
% * iccp_plot_crosscorr_pd_hw_ccc_pos_neg(ccpairs, showstats)
% * iccp_plot_strf_fr_pli_crosscorr_pos_neg(ccpairs)
% * iccp_plot_strf_similarity_crosscorr_pos_neg(ccpairs)
% * iccp_plot_strf_cf_q_latency_crosscorr_pos_neg(ccpairs) 
% * iccp_plot_strf_bmf_crosscorr_pos_neg(ccpairs)
% * iccp_plot_nonlinearity_asi_si_crosscorr_pos_neg(ccpairs)
% * iccp_plot_fiofit_crosscorr_pos_neg(ccpairs)


%% Saving data to a csv file
% * iccp_struct2csv(data,csvfile): data is from synchrony analysis, csvfile is user specified


%% Getting unique neurons and their values from ccpairs
% * [id,table] = iccp_ccpairs_to_table(ccpairs)






