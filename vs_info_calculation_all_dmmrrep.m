function vs_info_calculation_all_dmmrrep
all_sites = {'site19c','site19b','site18',...
    'site13',...
    'site12','site11','site4','site1'};

al = length(all_sites);
pathname = 'C:\Users\Victor\Documents\Research\20040217';

for i = 1:al
    newPathName = [pathname '\' all_sites{i}];
    cd(newPathName);
    d0 = dir('*-dmrrep-fs18115-final-resp.mat');
    load(d0.name);
    
    [info_nmse,sigma_nmse,nmse] = vs_info_calculation(resp,1,1);
    
end
end