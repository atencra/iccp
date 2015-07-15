function [position, cf, bw, q, latency] = get_iccpairs_strf_ptparams_folder
% get_iccpairs_strf_ptparams_folder BF, Q, latency from STRF files
% 
%    [position, bf, q, latency] = get_iccpairs_strf_ptparams_folder
% 
%    Reads through the current directory and finds files ending in 
%    '*-pairs-ptparams.mat'. These files hold the pure tone strf parameters
%    for cortical neurons. Next, the bf, spectral tuning (q), and
%    latency for pairs of neurons are extracted, and values 
%    for pairs of neurons are returned in position, bf, q, and latency, where:
% 
%    position : cortical depth of the pair
%    bf : best frequency
%    bw : bandwidth
%    q : spectral tuning
%    latency : in ms



cf = [];
bw = [];
q = [];
latency = [];
position = [];


d = dir('*-strfcmb-pairs-ptparams.mat');

for n = 1:length(d)

	filename = d(n).name;
	s = load(filename, 'ptparams');
   ptparams = s.ptparams;

   chan = [ptparams.chan];
   chan_unique = unique(chan);

   for i = 1:length(chan_unique)

      index = find( chan_unique(i) == chan );

      if ( length(index) > 1 )

         cmb = nchoosek(index, 2); % determine all possible pairs 

         [nr, nc] = size(cmb);

         for j = 1:nr

            index1 = cmb(j,1);
            index2 = cmb(j,2);

            chan3 = ptparams(index1).chan(1); 
            chan4 = ptparams(index2).chan(1); 
            model3 = ptparams(index1).model(1); 
            model4 = ptparams(index2).model(1); 

            pos = ptparams(index1).position;
            position = [position; pos];

            cf = [cf; ptparams(index1).cf ptparams(index2).cf];
            q = [q; ptparams(index1).q ptparams(index2).q];

            bw1 = abs( log2(ptparams(index1).fupper ./ ptparams(index1).flower  ) );
            bw2 = abs( log2(ptparams(index2).fupper ./ ptparams(index2).flower  ) );
            bw = [bw; bw1 bw2];

            latency1 = ptparams(index1).latency;
            latency2 = ptparams(index2).latency;

            latency = [latency; latency1 latency2];

         end % (for j)

      end % (if)

   end % (for i)

end % (for n)



return;







