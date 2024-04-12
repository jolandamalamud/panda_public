function SurrogateData = batchGenerateSurrogateData(Data,models,opt)
% 
% SurrogateData = batchModelFit(Data,modelsClass,options)
%
% Generate surrogate experimental data from fitted models contained in modelList
% to compare to original data. 
% 
% options.nSamples	how many samples per subject to generate (default: 100)
% 
% Quentin Huys, 2018 qhuys@cantab.net

if ~exist('opt'); opt = struct; end
if ~isfield(opt,'nSamples'); 		opt.nSamples = 100; 					end;
if ~isfield(opt,'resultsDir'); 	opt.resultsDir= [pwd filesep 'resultsDir'];end;

nModls = length(models);
Nsj = length(Data);

for mdl = 1:nModls

	% load EM-MAP parameters 
	try 
		R.(models(mdl).name) = load(sprintf('%s/%s.mat',opt.resultsDir,models(mdl).name));
		par = R.(models(mdl).name).E; 
	catch 
		fprintf('No fits for model %s found, not surrogate data generated\n',models(mdl).name);
		return 
	end

	fstr=str2func(models(mdl).name);	% turn variable into function call 
	doprior = 0; 							% add a prior for regularization or not
	mu = [];									% prior mean 
	nui = []; 								% prior variance 
	llopt.generatesurrogatedata=1;	% whether to generate surrogate data - here not 

	% generate structure for output 
	[foo,foo,foo2] = fstr(par(:,1),Data(1),mu,nui,doprior,llopt); 
	f = fieldnames(foo2)'; f{2,1} = {}; f = struct(f{:});
	dsurr0=repmat(f,1,opt.nSamples);

	fprintf('generating %i surrogate datasets from model %s\r',opt.nSamples,models(mdl).name);
	parfor sj=1:Nsj;
		dsurr=dsurr0;
		for ns=1:opt.nSamples
			[foo,foo,dsurr(ns)] = fstr(par(:,sj),Data(sj),mu,nui,doprior,llopt); 
		end
		SurrogateData(sj).(models(mdl).name) = dsurr; 
    end
    surr = vertcat(SurrogateData.(models(mdl).name));
    save([opt.resultsDir '/SurrogateData_' models(mdl).name '.mat'],'surr');
	fprintf('\n')
end

