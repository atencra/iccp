function dataFolder = iccp_fra_dir_from_spk(spk)
% iccp_fra_dir_from_spk Construct ICC relative directory structure from data struct
% 
%     dataFolder = iccp_fra_dir_from_spk(spk)
% 
%     Constructs the internal data path that is indicated by the fields in the 
%     struct spk. The experiment, site, depth, atten, and stim fields
%     are used to construct the path.
% 
%     All ICC data is stored in folders that look like:
% 
%     D:\InferiorColliculus\20040217\site14\2004-2-17-site14-5000um-30db-tc1
% 
%     This function constructs the 
%         \20040217\site14\2004-2-17-site14-5000um-30db-tc1
%     part using the data in spk.
% 
%     spk : struct array or struct element holding waveform and spike time 
%     data.


exp = spk(1).exp;
site = spk(1).site;
depth = spk(1).depth;
atten = spk(1).atten;
stim = spk(1).stim;


if ( strcmp(exp, '2004-2-17') )
    siteDir = '20040217';   
elseif ( strcmp(exp, '2003-11-24') )
    siteDir = '20031124';
elseif ( strcmp(exp, '2003-11-19') )
    siteDir = '20031119';
else
    error('No data for ICC on that date.');
end

dataFolder = fullfile(...
siteDir, ...
sprintf('site%.0f',site), ...
sprintf('%s-site%.0f-%.0fum-%.0fdb-%s',exp,site,depth,atten,stim) );


return;