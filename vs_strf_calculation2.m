function strf_calculation2
% strf_calculation2: Calculates STRF's from spk and trigger
%                   files, saves them in strf-new.mat files.

% Data taken from all sites with mdl_spk folders and dmrrep data, or sites
% 19c, 19b, 18, 13, 12, 11, 4, and 1.        

% Calls the function calculate_strf to calculate strfs.


envfile = 'C:\Users\Victor\Documents\Research\20040217\stimuli_20040217\dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min.spr';

all_sites = {'site19c','site19b','site18',...
    'site13',...
    'site12','site11','site4','site1'};

al = length(all_sites);
pathname = 'C:\Users\Victor\Documents\Research\20040217';

for i = 1:al
     newPathName = [pathname '\' all_sites{i} '\mdl_spk']
     cd(newPathName);
     d0 = dir('dmrrep');
     d1 = dir('dmr1');
     d2 = dir('dmr2');
     
     if ~isempty(d0) && ~isempty(d1)
         cd([newPathName '\dmr1']);
         [newPathName '\dmr1']
         d = dir('*-spk.mat');
         fileLoc = [cd '\' d.name];
         load(fileLoc,'spk');
         
         prefix = d.name;
         index = strfind(prefix,'spk');
         trigFileName = [prefix(1:index-1) 'trig.mat'];
         strfFileName = [prefix(1:index-1) '-new-strf.mat'];
         cd([pathname '\' all_sites{i}]);
         load(trigFileName,'trigger');
         
         
         strf = calculate_strf(spk, trigger, 0, envfile, 0, 0.1);
         save(strfFileName,'strf');
         
     end
     
     if ~isempty(d0) && ~isempty(d2)
         cd([newPathName '\dmr2']);
         [newPathName '\dmr2']
         d = dir('*-spk.mat');
         fileLoc = [cd '\' d.name];
         load(fileLoc,'spk');
         
         prefix = d.name;
         index = strfind(prefix,'spk');
         trigFileName = [prefix(1:index-1) 'trig.mat'];
         strfFileName = [prefix(1:index-1) '-new-strf.mat'];
         cd([pathname '\' all_sites{i}]);
         load(trigFileName,'trigger');
         
         strf = calculate_strf(spk, trigger, 0, envfile, 0, 0.1);
         save(strfFileName,'strf');
         
     end
     
     
     
end

         