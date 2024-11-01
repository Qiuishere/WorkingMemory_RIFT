function id = subject_id(prm, fields)
default('fields',{'number','code','age','sex'})
% input data
options     = struct('Resize','on','WindowStyle','modal');
id_input    = inputdlg(fields,prm.exp.name,1,...
                       repmat({''},numel(fields),1),options);
id          = cell2struct(id_input,fields);
% add 0 before 1-9
if isfield(id,'number') && str2double(id.('number'))<10                    %#> without dynamic structures: str2double(getfield(id,'number'))
    id.('number') = ['0' id.('number')];
end
% go to the exp prm.exp.path
% sep = '\';if ismac; sep = '/';end
cd(prm.exp.path);
   
% write the general xls
if ~ismac && isfield(id,'number')
xlswrite(prm.exp.name,fieldnames(id)' ,'participants','A1');  
xlswrite(prm.exp.name,struct2cell(id)','participants', ...
                  strcat('A',num2str(str2double(id.number)+1)));
else
disp('xlswrite not available for OSX or number not specified');
end

% % add other fields to id
% id.prm.exp.path   = prm.exp.path;
% id.prm.exp.name   = prm.exp.name;

% create the participant folder
content     = dir_ls;
found       = 0;
for k = 1:numel(content)
    if strcmp(content(k).name,'Results')
        found   = 1;
    end
end
if found==0
    mkdir('Results');   cd('Results');
    mkdir('Data');      cd('../'); 
end
cd('Results/Data');curdir = pwd;
dir(curdir)
if exist([pwd '\' id.number],'dir')>0 && ~strcmp(id.number,'99') % usually 99 is debug
    prompt = questdlg('this folder already exists, overwrite?', ...
                         'Overwrite', ...
                         'Yes','No','No');
    switch prompt
        case 'Yes'
            disp('overwriting...')
            rmdir([curdir '\' id.number],'s')
            cd(curdir)
        case 'No'
            error('check again your subjects''list and id.')
    end
end   
mkdir(id.number);
cd(id.number);
id.path     = pwd;
