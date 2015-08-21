function [cf_total, bw_total, q_total, latency_total] = iccp_combine_strf_ptparams
% iccp_combine_strf_ptparams Combine STRF pure tone data
%
% [cf_total, bw_total, q_total, latency_total] = iccp_combine_strf_ptparams
% --------------------------------------------------------------------------------
%
% Finds the *-strf-ptparams.mat files in a folder and combines all the 
% data. The files will have 4 variables: cf, bw, q, and latency.
%
% The variables will have been obtained fro ptparams struct arrays that are
% saved in *-strfcmb-pairs-ptparams.mat and have been processed by running
% the function 
%     [cf_total, bw_total, q_total, latency_total] = ...
%        iccp_combine_strf_ptparams;
%
% This is a helper that sets up a data base for later processing and
% plotting.
%
% caa 8/15/14


narginchk(0,0);

pfiles = dir('*-strf-ptparams.mat');
cf_total = [];
bw_total = [];
q_total = [];
latency_total = [];

for j = 1:length(pfiles)
   infile = pfiles(j).name;
   load(infile,'cf','bw','q','latency');
   cf_total = [cf_total; cf];
   clear cf;
   bw_total = [bw_total; bw];
   clear bw;
   q_total = [q_total; q];
   clear q;
   latency_total = [latency_total; latency];
   clear latency;
end % (for j)

return;




