function vs_plot_strf_pairs_multiple(strf,trigger,flim,tlim)
%plot_strf_pairs_multiple  Values from same channel cross-corr functions
%
% plot_strf_pairs_multiple
% Performs pairwise analysis on STRFs, then uses plot_strf_single to plot
% the paired STRFs in a 5x4 grid. The first and second columns are paired,
% and the third and fourth columns are paired.

if ( nargin == 0)
   d = dir('*-strfcmb-pairs.mat');
else
   d = 1;
   exp = strf(1).exp;
   site = strf(1).site;
   depth = strf(1).depth;
   ptitle = sprintf('%s-site%.0f-%.0fum',exp,site,depth);
end

for n = 1:length(d)
       
   if ( nargin == 0 )
       filename = d(n).name
       s = load(filename, 'strf', 'trigger');
       strf = s.strf;
       trigger = s.trigger;
   elseif ( nargin == 1 )
      error('You need to input the trigger argument.');
   elseif (nargin == 2)
      flim = [];
      tlim = [];
   elseif ( nargin == 3 )
      tlim = [];
   end
    
    chan = [strf.chan];
    chan_unique = unique(chan);
    
    strf_1 = [];
    strf_2 = [];
    
%     Finding Pairs
    for i = 1:length(chan_unique)
        
        index = find( chan_unique(i) == chan );
        
        if length(index)>1
            
            cmb = nchoosek(index, 2); % determine all possible pairwise combinations
            
            [nr, nc] = size(cmb);
            
            for j = 1:nr
                
                index1 = cmb(j,1);
                index2 = cmb(j,2);
                
                strf_1 = [strf_1; strf(index1)];
                strf_2 = [strf_2; strf(index2)];
                
            end % for j
            
        end % if
        
    end % for i
    
%     Plotting Pairs Next to Each Other
    
    L = length(strf_1);
    
    for m = 1:floor(L/10)
        
        figure;
        
        for k = 1:10
  
            if 2*k <= L 
                subplot(5,4,2*k-1)
                if ( nargin == 0 )
                   plot_strf_single(strf_1(5*(m-1)+k),trigger, [], [0 50])
                   box off;
                else
                   plot_strf_single(strf_1(5*(m-1)+k),trigger, flim, tlim)
                   box off;
                end
                
                subplot(5,4,2*k)
                if ( nargin == 0 )
                   plot_strf_single(strf_2(5*(m-1)+k),trigger, [], [0 50])
                   box off;
               else
                   plot_strf_single(strf_2(5*(m-1)+k),trigger, flim, tlim)
                   box off;
               end
            else
                subplot(5,4,2*k-1)
                subplot(5,4,2*k)
                
            end % if
           
        end % for k
        get(gcf,'position')

        set(gcf,'position',[100 100 665 661]);
        print_mfilename(mfilename);
        suptitle(ptitle);
        pause;
%         close;
    end % for m
    
end %for n

end