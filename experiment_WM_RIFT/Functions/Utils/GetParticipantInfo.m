function info = GetParticipantInfo(varargin)

% returns a structure with field of age, gender and hand.
    % Check if the input argument is not one of the defined attributes
%     if ~ismember(varargin{1}, attributes)
%         error('Invalid input argument. Only ''age'', ''gender'', or ''hand'' are allowed.');
%     end
    
    % Define the participant information based on the input argument
        title       =   sprintf('Information of Participant');
        for i = 1:numel(varargin)
            
            switch varargin{i}
                case 'age'
                    ques = {'Please enter your age'};
                    input      =   inputdlg( ques, title, [1 60]);
                    info.age = str2num(input{1});
                case 'gender'
                    options = {'female','male', 'other', 'prefer not to say'};
                    ques = strcat('Please enter your gender (1: ', options{1}, '; 2: ', options{2}, '; 3: ', options{3}, '; 4: ', options{4}, ')');
                    input      =   inputdlg( ques, title, [1 60]);
                    info.gender = options{str2num(input{1})};
                case 'hand'
                    options = {'Left','Right'};
                    ques = strcat('Which hand do you use the most often? (1: ', options{1}, '; 2: ', options{2}, ')');
                    input      =   inputdlg( ques, title, [1 60]);
                    info.hand = options{str2num(input{1})};
            end
                        
        end

end

