function Data=generateExampleDataset(Nsj,resultsDir, model, setparams, T)
% 
% Data = generateExampleDataset(Nsj)
% 
% Generate example dataset containing Nsj subjects for affective Go/Nogo task using the
% standard model llbaepqx.m
% 
% Quentin Huys 2018 www.quentinhuys.com 


fprintf('Generating example dataset for affective Go/Nogo task\n')

options.generatesurrogatedata=1; 

if isempty(T)
    T = 160;
end

fstr = str2func(model.name);

for sj=1:Nsj
	Data(sj).ID = sprintf('Subj %i',sj);
	
	Data(sj).a = zeros(1,T);					% preallocate space
	Data(sj).r = zeros(1,T);					% preallocate space
	rs = randperm(T);							% randomise stimuli 
	s = [1:4]'*ones(1,T/4);
	Data(sj).s = s(rs);							

	Data(sj).Nch = T; 							% length 
    
    if ~isempty(setparams)
        Data(sj).trueParam = setparams(:,sj);%mvnrnd(setparams.mu, setparams.sigma);
    else
        Data(sj).trueParam = randn(model.npar,1);
    end

	% generate choices A, state transitions S and rewards R 
	[foo,foo,dsurr] = fstr(Data(sj).trueParam,Data(sj),0,0,0,options); 
	Data(sj).a = dsurr.a;
	Data(sj).r = dsurr.r;
	Data(sj).trueModel = model.name;

end

% fprintf('Saved example dataset as Data.mat\n');
% save([resultsDir filesep 'Data.mat'],'Data');
end
