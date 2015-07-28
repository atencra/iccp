function nonlinearitypart3

% nonlinearitypart3: Generates images of STAS. Loads in the correct
% stimulus matrix and stimulus/locator data (contained in *-locator.mat),
% retrieves the STA from it with "get_sta_from_locator", and plots it using
% "imagesc".



load('dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8-matrix');
d = dir('*-locator.mat');
L = length(d)

for n = 1: L
   filename = d(n).name
   load(filename);
   
   for i = 1: length(spkloc)
       locator(:,i) = spkloc(i).locator;
       sta = get_sta_from_locator(locator(:,i),stimulus);
       if i == 1
            imagesc(sta);
            title(sprintf('%s, STA %.0f',filename, i));
            pause;
            close;
       end
   end
   
   figure;
   title(sprintf('Divider %.0f',n));

end