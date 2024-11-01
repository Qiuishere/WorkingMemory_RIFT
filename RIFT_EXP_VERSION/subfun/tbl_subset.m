function tbl = tbl_subset(tbl,varargin)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Extract rows from a dataset table
%                                                         D. Pascucci, EPFL
% Last update: 01.08.2021
%--------------------------------------------------------------------------
% INPUT
% - tbl:  dataset in table format
% - additional input pairs (e.g., tbl_subset(tbl,'condition','A'))
%   . each pair specify one colunm of the table and one value for the
%     selector, logical operator, nan selectors and functions can be used
%     as strings (see from line 25)
%--------------------------------------------------------------------------
% OUTPUT
% - tbl: table containing only the subset selected
%==========================================================================
if isstruct(varargin{1})
    varargin  = namedargs2cell(varargin{:});
end
names         = varargin(1:2:end);
values        = varargin(2:2:end);
assert(nnz(ismember(tbl.Properties.VariableNames,names))==numel(unique(names)), ...
    'some of the variables are not in the table');
i             = 1;
index         = zeros(size(tbl,1),1);
while i<=numel(varargin)/2
    assert(ischar(names{i}),'variable names should be char')
    v         = values{i};
    if ~isnumeric(v) && ~islogical(v) ...
            && nnz(ismember(v,'<>=~'))>0 && ~contains(v,'isnan')
        if nnz(ismember(v,'|'))>0
            v = strrep(v,'=','');
            temp = strsplit(v,'|');
            conditionText = '';
            for tt = 1:length(temp)
                switch tt
                    case length(temp)
                        conditionText =  strcat(conditionText,strcat(' x==', temp(tt)));
                    otherwise
                        conditionText =  strcat(conditionText,strcat(' x==', temp(tt),' | '));
                end
            end
            fun   = str2func(['@(x) ' conditionText{1} ]);
            index = index+fun(tbl.(names{i}));
        else
            fun   = str2func(['@(x) x' v]);
            index = index+fun(tbl.(names{i}));
        end
    elseif ~isnumeric(v) && contains(v,'isnan')
        fun   = str2func(['@(x) ' v '(x)' ]);
        index = index+fun(tbl.(names{i}));
    else        
        index = index+ismember(tbl.(names{i}),v);
    end
    i         = i+1;
end
% ensure all conditions are met
tbl           = tbl(mean(index,2)==(i-1),:);
