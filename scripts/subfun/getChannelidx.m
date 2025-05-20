function channelidx = getChannelidx(channelsLabel, channelList)
% Input:
% channelsLabel (cell) contain all the channels' label
% channleList (cell) contain all the channel to be selected
%
% Output:
% channelidx (array) corresponding number of the wanted channels 
    n             = numel(channelList);
    channelidx    = zeros(1,n);
    for i  = 1:n
        channelidx(i)  = find( strcmp( channelsLabel,channelList{i}));
    end
end