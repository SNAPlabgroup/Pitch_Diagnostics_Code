%% Set visit metadata
clear visit; 

prompts = {'Visit Date', 'Researcher', 'Subject ID', 'Age (yrs)', 'Experiment Group', 'Gender (M/F/X)'}; 
title = 'Visit Metadata'; 
dims = [1 20];
defaults = {date, 'Hauser', 'S', '', '', ''}; 

answer = inputdlg(prompts, title, dims, defaults); 

for i = length(answer)
    if isempty(answer{i})
        answer{i} = 'unknown'; 
    end
end

visit.room = 'LYLE3063'; 
visit.univ = 'Purdue'; 
visit.date = answer{1}; 
visit.researcher = answer{2};
visit.subj.ID = answer{3};
visit.subj.age = answer{4};
visit.subj.exp_group = answer{5};
visit.subj.gender = answer{6};

save('current_visit.mat', 'visit')