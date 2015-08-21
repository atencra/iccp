function iccp_process_fra_waveforms(spk,fra,trigger, chan,model,a,f)
% iccp_process_fra_waveforms Helper function to plot fra and spike shapes
% 
%     iccp_process_fra_waveforms(spk,fra,chan,model)
% 
%     spk : struct array of spike times, waveforms
% 
%     fra : struct array for tuning curves
% 
%     chan : integer, channel number to look at
% 
%     model : 1x2 vector, holds model numbers to look at
% 
%     This function collects some commands to make processing data go
%     faster. There is nothing sophisticated.
% 
%     One assumption: the path to the ICC data is assumed to be in:
%     'D:\InferiorColliculus';

narginchk(5,7);

if ( nargin == 5 )
    a = [];
    f = [];
end

if ( nargin == 6 )
    f = [];
end


p1 = 'D:\InferiorColliculus';

p2 = iccp_dir_from_spk(spk);

p3 = fullfile(p1,p2);

[e1] = ss_element_from_spk(spk, chan, model(1));
[e2] = ss_element_from_spk(spk, chan, model(2));

fra_plot_pairs_spk_fra(spk([e1 e2]),fra([e1 e2]));

ss_plot_waveform_trace_tone(p3, spk, [e1 e2], trigger,1627,fra,a,f);

return;