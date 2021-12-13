function PredSim(S,slv,pp,batchQueue)

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);

S = GetDefaultSettings(S);

% construct OpenSim model file name
OsimFileName = [S.subject '_' S.Foot.Model];
if strcmp(S.Foot.Scaling,'default')
    OsimFileName = [OsimFileName '_sd'];
elseif strcmp(S.Foot.Scaling,'custom')
    OsimFileName = [OsimFileName '_sc'];
elseif strcmp(S.Foot.Scaling,'personalised')
    OsimFileName = [OsimFileName '_sp'];
end
S.OsimFileName = OsimFileName;

% construct external function file name
ExternalFunc = OsimFileName;
if S.Foot.contactStiffnessFactor == 10
    ExternalFunc = [ExternalFunc '_cspx10'];
end
if S.Foot.contactSphereOffsetY
    ExternalFunc = [ExternalFunc '_oy'];
end
if S.Foot.contactSphereOffset1X
    ExternalFunc = [ExternalFunc '_o1x'];
end
S.ExternalFunc = ExternalFunc;



% build standardised names
[savename, casfuncfol] = getSavename(S);
S.CasadiFunc_Folders = casfuncfol;
S.savename = savename;


if batchQueue
    if (exist([pathRepo '/Results/batchQ.mat'],'file')==2) 
        load([pathRepo '/Results/batchQ.mat'],'batchQ');
    else
        batchQ.(S.savename) = struct('S',[]);
    end
    batchQ.(S.savename).S = S;
    
else
    % make folder to store results if it doesn't exist
    pathResults = fullfile([pathRepo '/Results'],S.ResultsFolder);
    if ~isfolder(pathResults)
        mkdir(pathResults);
    end

    % Create the casadifunctions if they do not exist yet
    if ~isfolder([pathRepo '\CasADiFunctions\' S.CasadiFunc_Folders])
        disp('Creating casadifunctions...');
        CreateCasadiFunctions(pathRepo,S);
        disp('...casadifunctions created');
    end

end

%% Run

if batchQueue
    batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_FootModel';
    batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_FootModel';
end
if slv        % run the optimization
    f_PredSim_Gait92_FootModel(S);
end
if pp           % post-proces simulation results
    f_LoadSim_Gait92_FootModel(S.ResultsFolder,S.savename);
end


if batchQueue
    save([pathRepo '/Results/batchQ.mat'],'batchQ');
end













