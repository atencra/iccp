function [iccdata, a1data] = plot_pairs_icc_a1_spktype_strf(a1data, iccdata)
% plot_pairs_strf_cf_q_latency - compare pure tone params of pairs
%
% 
% plot_pairs_strf_cf_q_latency
% --------------------------------------------------------
% The function finds the neurons that were recorded from the 
% same channel, and then compares the pure tone parameters
% for these neurons.
%
% The analysis tells us the variability of the parameters in
% a locally confined area of cortex.
%
% caa 4/1/11
%
% plot_pairs_strf_cf_q_latency;


if ( nargin == 0 )
   iccdata = [];
   a1data = [];

   [position1, cf, q, latency] = get_pairs_spktype_strf_ptparams_folder;
   [position2, fr, rpi, spi] = get_pairs_spktype_strf_params_folder;
   [position3, si] = get_pairs_spktype_strf_similarity_folder;
   [position4, fiosi, fioasi] = get_pairs_spktype_fio_nonlinearity_asi_si_folder;

   [size(cf) size(q) size(latency)] 
   [size(fr) size(rpi) size(spi)]
   [size(si)]
   [size(fiosi) size(fioasi)]
   [size(position1) size(position2) size(position3) size(position4)]

   a1data.position1 = position1;
   a1data.cf = cf;
   a1data.q = q;
   a1data.latency = latency;

   a1data.position2 = position2;
   a1data.fr = fr;
   a1data.rpi = rpi;
   a1data.spi = spi;

   a1data.position3 = position3;
   a1data.si = si;

   a1data.position4 = position4;
   a1data.fiosi = fiosi;
   a1data.fioasi = fioasi;



   position = [ccpairs.position];
   position = position(:);

   cf1 = [ccpairs.cf1];
   cf2 = [ccpairs.cf2];
   cf = [cf1(:) cf2(:)];

   q1 = [ccpairs.q1];
   q2 = [ccpairs.q2];
   q = [q1(:) q2(:)];

   latency1 = [ccpairs.latency1];
   latency2 = [ccpairs.latency2];
   latency = [latency1(:) latency2(:)];

   fr1 = [ccpairs.fr1];
   fr2 = [ccpairs.fr2];
   fr = [fr1(:) fr2(:)];

   rpi1 = [ccpairs.pli1];
   rpi2 = [ccpairs.pli2];
   rpi = [rpi1(:) rpi2(:)];

   spi1 = [ccpairs.sepindex1];
   spi2 = [ccpairs.sepindex2];
   spi = [spi1(:) spi2(:)];

   si = [ccpairs.strfsimilarity];
   si = si(:);

   asi1 = [ccpairs.fioasi1];
   asi2 = [ccpairs.fioasi2];
   fioasi = [asi1(:) asi2(:)];
   fiosi = [ccpairs.fiosimilarity];
   fiosi = fiosi(:);

   iccdata.position = position(:);
   iccdata.cf = cf;
   iccdata.q = q;
   iccdata.latency = latency;
   iccdata.fr = fr;
   iccdata.rpi = rpi;
   iccdata.spi = spi;
   iccdata.si = si;
   iccdata.fiosi = fiosi;
   iccdata.fioasi = fioasi;


end;


close all;

% plot_pairs_cf_q_latency_granular_diff_bar(a1data, iccdata);

% plot_pairs_fr_pli_spi_granular_bar(a1data, iccdata);

plot_pairs_strf_similarity_granular_bar(a1data, iccdata);

% plot_pairs_fio_similarity_granular_bar(a1data, iccdata);

% plot_pairs_fio_asymmetry_granular_bar(a1data, iccdata);

return;




function plot_pairs_cf_q_latency_granular_diff_bar(a1data, iccdata) 
% plot_pairs_population_cf_q_latency - neuron pairs mtf parameters across layer
%
% plot_pairs_population_mtf(btmf, bsmf, twc3db, swc3db)
% -----------------------------------------------------------------------
%
%
%
% caa 4/21/11

position = a1data.position1;
indexSupra = find(position < 600);
indexGran = find(position > 700 & position < 1100);
indexInfra = find(position > 1200);

% cf
datastr(1).position = position;
datastr(1).indexSupra = indexSupra;
datastr(1).indexGran = indexGran;
datastr(1).indexInfra = indexInfra;
datastr(1).data = a1data.cf;
datastr(1).iccdata = iccdata.cf;
datastr(1).datatype = 'cf';
datastr(1).xmin = 5;
datastr(1).xmax = 30;
datastr(1).ymin = 5;
datastr(1).ymax = 30;
datastr(1).xscale = 'log';
datastr(1).yscale = 'log';
datastr(1).edges = linspace(0, 0.5, 11);
datastr(1).tick = [5 10 20 30];
datastr(1).ticklabel = [5 10 20 30];
datastr(1).xlabel1 = 'BF (kHz)';
datastr(1).xlabel2 = 'Mean BF Diff (oct)';

% q
datastr(2).position = position;
datastr(2).indexSupra = indexSupra;
datastr(2).indexGran = indexGran;
datastr(2).indexInfra = indexInfra;
datastr(2).data = a1data.q;
datastr(2).iccdata = iccdata.q;
datastr(2).datatype = 'q';
datastr(2).xmin = 0.25;
datastr(2).xmax = 10;
datastr(2).ymin = 0.25;
datastr(2).ymax = 10;
datastr(2).xscale = 'log';
datastr(2).yscale = 'log';
datastr(2).edges = linspace(0, 5, 11);
datastr(2).tick = [0.25 0.5 1 2 4 8];
datastr(2).ticklabel = [0.25 0.5 1 2 4 8];
datastr(2).xlabel1 = 'Q';
datastr(2).xlabel2 = 'Mean Q Diff (oct)';

% latency
datastr(3).position = position;
datastr(3).indexSupra = indexSupra;
datastr(3).indexGran = indexGran;
datastr(3).indexInfra = indexInfra;
datastr(3).data = a1data.latency;
datastr(3).iccdata = iccdata.latency;
datastr(3).datatype = 'latency';
datastr(3).xmin = 0;
datastr(3).xmax = 30;
datastr(3).ymin = 0;
datastr(3).ymax = 30;
datastr(3).xscale = 'linear';
datastr(3).yscale = 'linear';
datastr(3).edges = linspace(0, 30, 11);
datastr(3).tick = [0:10:50];
datastr(3).ticklabel = [0:10:50];
datastr(3).xlabel1 = 'Latency (ms)';
datastr(3).xlabel2 = 'Mean Latency Diff (ms)';


ymaxdiff = [1 1 1];
ymaxdiff2 = [0.25 0.8 2];

for i = 1:length(datastr)

   data = datastr(i).data;
   idata = datastr(i).iccdata;
   indSupra = datastr(i).indexSupra;
   indGran = datastr(i).indexGran;
   indInfra = datastr(i).indexInfra;
   datatype = datastr(i).datatype;
   xmin = datastr(i).xmin;
   xmax = datastr(i).xmax;
   ymin = datastr(i).ymin;
   ymax = datastr(i).ymax;
   xscale = datastr(i).xscale;
   yscale = datastr(i).yscale;
   edges = datastr(i).edges;
   tick = datastr(i).tick;
   ticklabel = datastr(i).ticklabel;
   xlabel1 = datastr(i).xlabel1;
   xlabel2 = datastr(i).xlabel2;


   if ( strcmp(datatype, 'cf') )
      datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
      idatadiff = abs( log2( idata(:,1) ./ idata(:,2) ) );
   elseif ( strcmp(datatype, 'latency') )
      datadiff = abs( data(:,1) - data(:,2) );
      idatadiff = abs( idata(:,1) - idata(:,2) );
   else % it must be 'q'
      datadiff = abs( log2( data(:,1) ./ data(:,2) ) );
      idatadiff = abs( log2( idata(:,1) ./ idata(:,2) ) );
   end

   dxSupra = datadiff(indSupra);
   dxGran = datadiff(indGran);
   dxInfra = datadiff(indInfra);

   Pop = simple_stats(datadiff);
   Supra = simple_stats(dxSupra);
   Gran = simple_stats(dxGran);
   Infra = simple_stats(dxInfra);

   iPop = simple_stats(idatadiff);

   [x_cdf_pop, y_cdf_pop] = cumdist(datadiff);
   [x_cdf_supra, y_cdf_supra] = cumdist(dxSupra);
   [x_cdf_gran, y_cdf_gran] = cumdist(dxGran);
   [x_cdf_icc, y_cdf_icc] = cumdist(idatadiff);


   figure;
   subplot(2,4,[1 6]);
   hold on;
   cmap = brewmaps('orrd', 4);
   cmap = cmap(1:3,:);
%    plot(x_cdf_gran, y_cdf_gran, '-', 'color', cmap(2,:), 'linewidth', 2);
%    plot(x_cdf_icc, y_cdf_icc, '-', 'color', cmap(3,:), 'linewidth', 2);
   plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
   plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
   ytick = linspace(0,1,5);
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   tickpref;
   ylabel('Proportion');
   xlabel(xlabel2);
   legend('Gran', 'ICC');
%    xlim([0 1]);
   ylim([0 1.05]);
   [h, pval] = kstest2(datadiff, idatadiff);
   title(sprintf('p = %.4f', pval));


   subplot(2,4,[4 8]);
   hold on;
   v = [0 Gran.mn iPop.mn 0];
   hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
   children = get(hb,'children');
   cdata = repmat(1:numel(v),5,1);
   cdata = [cdata(:);nan];
   set(children,'facevertexcdata',cdata);
   cmap = brewmaps('orrd', 4);
   cmap = cmap(1:3,:);
   colormap([0.6 0.6 0.6; 0 0 0]);
   plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
   plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
   plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
   plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
   box off;
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   xtick = 1:2;
   xticklabel = {'Gran', 'ICC'};
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   xlabel('Station');
   ytick = linspace(0,ymaxdiff2(i),3);
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   xlim([0 3]);
   ylim([0 ymaxdiff2(i)]);
   ylabel(xlabel2);
   set(gcf,'position', [628   523   438   156]);
   print_mfilename(mfilename);

end % for i


return;



function plot_pairs_fr_pli_spi_granular_bar(a1data, iccdata)

position = a1data.position1;
indexSupra = find(position < 600);
indexGran = find(position > 700 & position < 1100);
indexInfra = find(position > 1200);


% firing rate
datatype(1).position = position;
datatype(1).indexSupra = indexSupra;
datatype(1).indexGran = indexGran;
datatype(1).indexInfra = indexInfra;
datatype(1).data = a1data.fr;
datatype(1).iccdata = iccdata.fr;
datatype(1).xmin = 0.1;
datatype(1).xmax = 30;
datatype(1).ymin = 0.1;
datatype(1).ymax = 30;
datatype(1).edges = 0:2.5:20;
datatype(1).tick = [0.1 1 10 30];
datatype(1).xscale = 'log';
datatype(1).yscale = 'log';
datatype(1).xlabel1 = 'Firing Rate';
datatype(1).xlabel2 = 'Firing Rate Difference (Hz)';



% phase-locking index
datatype(2).position = position;
datatype(2).indexSupra = indexSupra;
datatype(2).indexGran = indexGran;
datatype(2).indexInfra = indexInfra;
datatype(2).data = a1data.rpi;
datatype(2).iccdata = iccdata.rpi;
datatype(2).xmin = 0.01;
datatype(2).xmax = 1;
datatype(2).ymin = 0.01;
datatype(2).ymax = 1;
datatype(2).edges = 0:0.05:0.4;
datatype(2).tick = [0.01 0.1 1];
datatype(2).xscale = 'log';
datatype(2).yscale = 'log';
datatype(2).xlabel1 = 'Response Precision Index';
datatype(2).xlabel2 = 'RPI Difference';


% separability index
datatype(3).position = position;
datatype(3).indexSupra = indexSupra;
datatype(3).indexGran = indexGran;
datatype(3).indexInfra = indexInfra;
datatype(3).data = a1data.spi;
datatype(3).iccdata = iccdata.spi;
datatype(3).xmin = 0.2;
datatype(3).xmax = 1;
datatype(3).ymin = 0.2;
datatype(3).ymax = 1;
datatype(3).edges = 0:0.075:0.6;
datatype(3).tick = [0:0.2:1];
datatype(3).xscale = 'log';
datatype(3).yscale = 'log';
datatype(3).xlabel1 = 'Separability Index';
datatype(3).xlabel2 = 'SI Difference';

ymaxdiff = [0.53 0.7 0.6];
ytickmaxdiff = [0.45 0.6 3];

ytickmaxdiff2 = [6 0.06 3];
ymaxdiff2 = [8 0.07 3];

xlimcdf = [100 0.5 0.6];

cmap = cschemes('blues', 4);
cmap = cmap(1:3,:);

cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);

marker = {'o', 's', '^'};

for i = 1:2 %length(datatype)

   data = datatype(i).data;
   idata = datatype(i).iccdata;
   datadiff = abs( data(:,1) - data(:,2) );
   idatadiff = abs( idata(:,1) - idata(:,2) );
   indSupra = datatype(i).indexSupra;
   indGran = datatype(i).indexGran;
   indInfra = datatype(i).indexInfra;
   xmin = datatype(i).xmin;
   xmax = datatype(i).xmax;
   ymin = datatype(i).ymin;
   ymax = datatype(i).ymax;
   xscale = datatype(i).xscale;
   yscale = datatype(i).yscale;
   edges = datatype(i).edges;
   tick = datatype(i).tick;
   xlabel1 = datatype(i).xlabel1;
   xlabel2 = datatype(i).xlabel2;

   dxSupra = datadiff(indSupra);
   dxGran = datadiff(indGran);
   dxInfra = datadiff(indInfra);

   Pop = simple_stats(datadiff);
   Supra = simple_stats(dxSupra);
   Gran = simple_stats(dxGran);
   Infra = simple_stats(dxInfra);
   iPop = simple_stats(idatadiff);

   [x_cdf_pop, y_cdf_pop] = cumdist(datadiff);
   [x_cdf_supra, y_cdf_supra] = cumdist(dxSupra);
   [x_cdf_gran, y_cdf_gran] = cumdist(dxGran);
   [x_cdf_icc, y_cdf_icc] = cumdist(idatadiff);


   figure;
   subplot(2,4,[1 6]);
   hold on;
   cmap = brewmaps('orrd', 4);
   cmap = cmap(1:3,:);
%    plot(x_cdf_gran, y_cdf_gran, '-', 'color', cmap(2,:), 'linewidth', 2);
%    plot(x_cdf_icc, y_cdf_icc, '-', 'color', cmap(3,:), 'linewidth', 2);
   plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
   plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
   ytick = linspace(0,1,5);
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   tickpref;
   ylabel('Proportion');
   xlabel(xlabel2);
   legend('Gran', 'ICC');
   xlim([0 xlimcdf(i)]);
   ylim([0 1.05]);
   [h, pval] = kstest2(datadiff, idatadiff);
   title(sprintf('p = %.4f', pval));

   subplot(2,4,[4 8]);
   hold on;
   v = [0 Gran.mn iPop.mn 0];
   hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
   children = get(hb,'children');
   cdata = repmat(1:numel(v),5,1);
   cdata = [cdata(:);nan];
   set(children,'facevertexcdata',cdata);
   cmap = brewmaps('orrd', 4);
   cmap = cmap(1:3,:);
   colormap([0.6 0.6 0.6; 0 0 0]);
   plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
   plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
   plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
   plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
   box off;
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   xtick = 1:2;
   xticklabel = {'Gran', 'ICC'};
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   xlabel('Station');
   ytick = linspace(0,ymaxdiff2(i),3);
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   xlim([0 3]);
   ylim([0 ymaxdiff2(i)]);
   ylabel(xlabel2);
   set(gcf,'position', [628   523   438   156]);
   print_mfilename(mfilename);

end % for i



return;



function plot_pairs_strf_similarity_granular_bar(a1data, iccdata)

position = a1data.position3;
indexSupra = find(position <= 600);
indexGran = find(position >= 700 & position <= 1100);
indexInfra = find(position >= 1200);

si = a1data.si;
siSupra = si(indexSupra);
siGran = si(indexGran);
siInfra = si(indexInfra);
isi = iccdata.si;

Pop = simple_stats(si);
Supra = simple_stats(siSupra);
Gran = simple_stats(siGran);
Infra = simple_stats(siInfra);
iPop = simple_stats(isi);

[x_cdf_pop, y_cdf_pop] = cumdist(si);
[x_cdf_supra, y_cdf_supra] = cumdist(siSupra);
[x_cdf_gran, y_cdf_gran] = cumdist(siGran);
[x_cdf_icc, y_cdf_icc] = cumdist(isi);



figure;
subplot(2,4,[1 6]);
hold on;
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
% plot(x_cdf_gran, y_cdf_gran, '-', 'color', cmap(2,:), 'linewidth', 2);
% plot(x_cdf_icc, y_cdf_icc, '-', 'color', cmap(3,:), 'linewidth', 2);
plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
ytick = linspace(0,1,5);
set(gca,'ytick', ytick, 'yticklabel', ytick);
tickpref;
ylabel('Proportion');
xlabel('STRF Similarity');
legend('Gran', 'ICC');
xlim([0 1]);
ylim([0 1.05]);
[h, pval] = kstest2(si, isi);
title(sprintf('p = %.4f', pval));

subplot(2,4,[4 8]);
hold on;
v = [0 Gran.mn iPop.mn 0];
hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
children = get(hb,'children');
cdata = repmat(1:numel(v),5,1);
cdata = [cdata(:);nan];
set(children,'facevertexcdata',cdata);
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
colormap([0.6 0.6 0.6; 0 0 0]);
plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xtick = 1:2;
xticklabel = {'Gran', 'ICC'};
set(gca,'xtick', xtick, 'xticklabel', xticklabel);
xlabel('Station');
ytick = linspace(0,0.6,3);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([0 3]);
ylim([0 0.6]);
ylabel('STRF Similarity');
set(gcf,'position', [628   523   438   156]);
print_mfilename(mfilename);










si = 1 - a1data.si;
siSupra = si(indexSupra);
siGran = si(indexGran);
siInfra = si(indexInfra);
isi = 1 - iccdata.si;

Pop = simple_stats(si);
Supra = simple_stats(siSupra);
Gran = simple_stats(siGran);
Infra = simple_stats(siInfra);
iPop = simple_stats(isi);

[x_cdf_pop, y_cdf_pop] = cumdist(si);
[x_cdf_supra, y_cdf_supra] = cumdist(siSupra);
[x_cdf_gran, y_cdf_gran] = cumdist(siGran);
[x_cdf_icc, y_cdf_icc] = cumdist(isi);



figure;
subplot(2,4,[1 6]);
hold on;
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
ytick = linspace(0,1,5);
set(gca,'ytick', ytick, 'yticklabel', ytick);
tickpref;
ylabel('Proportion');
xlabel('1 - STRF Similarity');
legend('Gran', 'ICC');
xlim([0 1]);
ylim([0 1.05]);
[h, pval] = kstest2(si, isi);
title(sprintf('p = %.4f', pval));


subplot(2,4,[4 8]);
hold on;
v = [0 Gran.mn iPop.mn 0];
hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
children = get(hb,'children');
cdata = repmat(1:numel(v),5,1);
cdata = [cdata(:);nan];
set(children,'facevertexcdata',cdata);
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
colormap([0.6 0.6 0.6; 0 0 0]);
plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xtick = 1:2;
xticklabel = {'Gran', 'ICC'};
set(gca,'xtick', xtick, 'xticklabel', xticklabel);
xlabel('Station');
ytick = linspace(0,0.1,3);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([0 3]);
ylim([0 1]);
ylabel('1 - STRF Similarity');
set(gcf,'position', [628   523   438   156]);
print_mfilename(mfilename);


return;



function plot_pairs_fio_similarity_granular_bar(a1data, iccdata)


position = a1data.position4;
indexSupra = find(position <= 600);
indexGran = find(position >= 700 & position <= 1100);
indexInfra = find(position >= 1200);

si = a1data.fiosi;
siSupra = si(indexSupra);
siGran = si(indexGran);
siInfra = si(indexInfra);
isi = iccdata.fiosi;

Pop = simple_stats(si);
Supra = simple_stats(siSupra);
Gran = simple_stats(siGran);
Infra = simple_stats(siInfra);
iPop = simple_stats(isi);

[x_cdf_pop, y_cdf_pop] = cumdist(si);
[x_cdf_supra, y_cdf_supra] = cumdist(siSupra);
[x_cdf_gran, y_cdf_gran] = cumdist(siGran);
[x_cdf_icc, y_cdf_icc] = cumdist(isi);



figure;
subplot(2,4,[1 6]);
hold on;
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
% plot(x_cdf_gran, y_cdf_gran, '-', 'color', cmap(2,:), 'linewidth', 2);
% plot(x_cdf_icc, y_cdf_icc, '-', 'color', cmap(3,:), 'linewidth', 2);
plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
ytick = linspace(0,1,5);
set(gca,'ytick', ytick, 'yticklabel', ytick);
tickpref;
ylabel('Proportion');
xlabel('Nonlinearity Similarity');
legend('Gran', 'ICC');
xlim([0.4 1]);
ylim([0 1.05]);
[h, pval] = kstest2(si, isi);
title(sprintf('p = %.4f', pval));


subplot(2,4,[4 8]);
hold on;
v = [0 Gran.mn iPop.mn 0];
hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
children = get(hb,'children');
cdata = repmat(1:numel(v),5,1);
cdata = [cdata(:);nan];
set(children,'facevertexcdata',cdata);
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
colormap([0.6 0.6 0.6; 0 0 0]);
plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xtick = 1:2;
xticklabel = {'Gran', 'ICC'};
set(gca,'xtick', xtick, 'xticklabel', xticklabel);
xlabel('Station');
ytick = linspace(0,1,3);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([0 3]);
ylim([0 1]);
ylabel('Nonlinearity Similarity');
set(gcf,'position', [628   523   438   156]);
print_mfilename(mfilename);









si = a1data.fiosi;
si = 1 - si;
siSupra = si(indexSupra);
siGran = si(indexGran);
siInfra = si(indexInfra);
isi = iccdata.fiosi;
isi = 1 - isi;

Pop = simple_stats(si);
Supra = simple_stats(siSupra);
Gran = simple_stats(siGran);
Infra = simple_stats(siInfra);
iPop = simple_stats(isi);

[x_cdf_pop, y_cdf_pop] = cumdist(si);
[x_cdf_supra, y_cdf_supra] = cumdist(siSupra);
[x_cdf_gran, y_cdf_gran] = cumdist(siGran);
[x_cdf_icc, y_cdf_icc] = cumdist(isi);



figure;
subplot(2,4,[1 6]);
hold on;
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
% plot(x_cdf_gran, y_cdf_gran, '-', 'color', cmap(2,:), 'linewidth', 2);
% plot(x_cdf_icc, y_cdf_icc, '-', 'color', cmap(3,:), 'linewidth', 2);
plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
ytick = linspace(0,1,5);
set(gca,'ytick', ytick, 'yticklabel', ytick);
tickpref;
ylabel('Proportion');
xlabel('1 - Nonlinearity Similarity');
legend('Gran', 'ICC');
xlim([0 0.6]);
ylim([0 1.05]);
[h, pval] = kstest2(si, isi);
title(sprintf('p = %.4f', pval));


subplot(2,4,[4 8]);
hold on;
v = [0 Gran.mn iPop.mn 0];
hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
children = get(hb,'children');
cdata = repmat(1:numel(v),5,1);
cdata = [cdata(:);nan];
set(children,'facevertexcdata',cdata);
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
colormap([0.6 0.6 0.6; 0 0 0]);
plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xtick = 1:2;
xticklabel = {'Gran', 'ICC'};
set(gca,'xtick', xtick, 'xticklabel', xticklabel);
xlabel('Station');
ytick = linspace(0,0.06,3);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([0 3]);
ylim([0 0.06]);
ylabel('1 - Nonlinearity Similarity');
set(gcf,'position', [628   523   438   156]);
print_mfilename(mfilename);

return;



function plot_pairs_fio_asymmetry_granular_bar(a1data, iccdata)

position = a1data.position4;
indexSupra = find(position < 600);
indexGran = find(position > 700 & position < 1100);
indexInfra = find(position > 1200);

asi = a1data.fioasi;
datadiff = abs( asi(:,1) - asi(:,2) );
dxSupra = datadiff(indexSupra);
dxGran = datadiff(indexGran);
dxInfra = datadiff(indexInfra);

iasi = iccdata.fioasi;
idatadiff = abs( iasi(:,1) - iasi(:,2) );

Pop = simple_stats(datadiff);
Supra = simple_stats(dxSupra);
Gran = simple_stats(dxGran);
Infra = simple_stats(dxInfra);
iPop = simple_stats(idatadiff);

[x_cdf_pop, y_cdf_pop] = cumdist(datadiff);
[x_cdf_supra, y_cdf_supra] = cumdist(dxSupra);
[x_cdf_gran, y_cdf_gran] = cumdist(dxGran);
[x_cdf_icc, y_cdf_icc] = cumdist(idatadiff);


figure;
subplot(2,4,[1 6]);
hold on;
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
% plot(x_cdf_gran, y_cdf_gran, '-', 'color', cmap(2,:), 'linewidth', 2);
% plot(x_cdf_icc, y_cdf_icc, '-', 'color', cmap(3,:), 'linewidth', 2);
plot(x_cdf_gran, y_cdf_gran, '-', 'color', 0.6*ones(1,3), 'linewidth', 2);
plot(x_cdf_icc, y_cdf_icc, '-', 'color', 'k', 'linewidth', 2);
ytick = linspace(0,1,5);
set(gca,'ytick', ytick, 'yticklabel', ytick);
tickpref;
ylabel('Proportion');
xlabel('ASI Difference');
legend('Gran', 'ICC');
xlim([0 0.3]);
ylim([0 1.05]);
[h, pval] = kstest2(datadiff, idatadiff);
title(sprintf('p = %.4f', pval));



subplot(2,4,[4 8]);
hold on;
v = [0 Gran.mn iPop.mn 0];
hb = bar([0:3], [0 Gran.mn iPop.mn 0]);
children = get(hb,'children');
cdata = repmat(1:numel(v),5,1);
cdata = [cdata(:);nan];
set(children,'facevertexcdata',cdata);
cmap = brewmaps('orrd', 4);
cmap = cmap(1:3,:);
colormap([0.6 0.6 0.6; 0 0 0]);
plot([1 1], [Gran.mn Gran.mn+Gran.se], 'k-');
plot([0.8 1.2], [Gran.mn+Gran.se Gran.mn+Gran.se], 'k-');
plot([2 2], [iPop.mn iPop.mn+iPop.se], 'k-');
plot([1.8 2.2], [iPop.mn+iPop.se iPop.mn+iPop.se], 'k-');
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xtick = 1:2;
xticklabel = {'Gran', 'ICC'};
set(gca,'xtick', xtick, 'xticklabel', xticklabel);
xlabel('Station');
ytick = linspace(0,0.075,3);
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlim([0 3]);
ylim([0 0.075]);
ylabel('Mean ASI Difference');
set(gcf,'position', [628   523   438   156]);
print_mfilename(mfilename);

return;









