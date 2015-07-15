function strf_calculation
all_sites = {'site19c','site19b','site19','site18c','site18b','site18',...
    'site17b','site17','site16b','site16','site15','site14','site13',...
    'site12','site11','site10b','site10','site9','site8','site7','site6',...
    'site5','site4','site3b','site3','site2','site1'};
internal_folders = {'mdl_spk_individual_sort','mdl_spk'};
stimuli = {'dmr1','dmr2','ct','tc1','dmrrep','dmr1dmr2dmrrep'}; 

pathname = cd;
al = length(all_sites);
inf = length(internal_folders);

for i = 1:al
    
     for j = 1:inf
        
        for k = 1:length(stimuli)
            
            if j ==1
                newPathName = [pathname '\' all_sites{i} '\' all_sites{i} '_' internal_folders{j} '\' stimuli{k}]
                cd(newPathName);
                
            else
                newPathName = [pathname '\' all_sites{i} '\' internal_folders{j} '\' stimuli{k}]
                cd(newPathName);
            end
            
            d = dir('*.spk');
            if ~isempty(d)
                [spk, outfile] = saveSpikeSort(newPathName, stimuli{k}, 150);
                save(outfile, 'spk');
            end
            
        end
        
     end
    
end
