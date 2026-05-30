clear
% Define the report file name
reportFile = 'False_positive_TIC.txt';

% Open the file for writing (will overwrite if it exists)
fid = fopen(reportFile, 'w');

fprintf(fid, 'To check the presence of false positive TICs\n');
fprintf(fid, 'Generated on: %s\n\n', datestr(now));

% For the AGORA DB TICs and models
path2AGORATICs = '..\AGORA_TIC_Results\';
fileNamesAGORA = dir(path2AGORATICs);
fileNamesAGORA = {fileNamesAGORA(3:end).name}';
path2AGORAmodels = 'D:\OneDrive_moved_files\Microbial_community_gap_filling\With_AGORA2\AGORA2\';

% For the BiGG DB TICs and models
path2BIGGTICs = '..\BiggBlkdResults\';
fileNamesBIGG = dir(path2BIGGTICs);
fileNamesBIGG = {fileNamesBIGG(3:end).name}';
path2BIGGmodels = 'D:\OneDrive - smail.iitm.ac.in\SprintCore\MinReact\BiggModelsMat\';
fileNamesBIGG = fileNamesBIGG(~contains(fileNamesBIGG,'Incomp'));

% Verifying the absence of false postive TICs in the AGORA models
N = numel(fileNamesAGORA); % Total number of models
for k = 1:N
    % loading the model
    load([path2AGORAmodels,fileNamesAGORA{k}]);
    % loading the TICs
    load([path2AGORATICs,fileNamesAGORA{k}]);
    
    for t =1:numel(TICs)
        ticmodel = removeRxns(model,setdiff(model.rxns,TICs{t}));
        a= fastcc(ticmodel,1e-4,0);
        % checking for the consistency
        if numel(ticmodel.rxns)~=numel(a)
            msg = sprintf('The model, %s has a TIC that is inconsistent', fileNamesAGORA{k});
            fprintf(fid, msg);
        end
        % checking for minimality-I
        if rank(full(ticmodel.S))~=numel(ticmodel.rxns)-1
            msg = sprintf('The model, %s has a TIC that is not minimal', fileNamesAGORA{k});
            fprintf(fid, msg);
        end
        
        % checking for minimality-II
        rand_rxn = randsample(numel(ticmodel.rxns),1);
        ticmodel.lb(rand_rxn)=0;ticmodel.ub(rand_rxn)=0;
        a= fastcc(ticmodel,1e-4,0);
        if ~isempty(a)
            msg = sprintf('The model, %s has a TIC that is not minimal', fileNamesAGORA{k});
            fprintf(fid, msg);
        end
    end
end

% Verifying the absence of false postive TICs in the BIGG models
N = numel(fileNamesBIGG); % Total number of models
for k = 1:N
    % loading the model
    load([path2BIGGmodels,fileNamesBIGG{k}]);
    % loading the TICs
    load([path2BIGGTICs,fileNamesBIGG{k}]);
    
    for t =1:numel(TICs)
        ticmodel = removeRxns(model,setdiff(model.rxns,TICs{t}));
        a= fastcc(ticmodel,1e-4,0);
        % checking for the consistency
        if numel(ticmodel.rxns)~=numel(a)
            msg = sprintf('The model, %s has a TIC that is inconsistent', fileNamesAGORA{k});
            fprintf(fid, msg);
        end
        % checking for minimality-I
        if rank(ticmodel.S)~=numel(ticmodel.rxns)-1
            msg = sprintf('The model, %s has a TIC that is not minimal', fileNamesAGORA{k});
            fprintf(fid, msg);
        end
        
        % checking for minimality-II
        rand_rxn = randsample(numel(ticmodel.rxns),1);
        ticmodel.lb(rand_rxn)=0;ticmodel.ub(rand_rxn)=0;
        a= fastcc(ticmodel,1e-4,0);
        if ~isempty(a)
            msg = sprintf('The model, %s has a TIC that is not minimal', fileNamesAGORA{k});
            fprintf(fid, msg);
        end
    end
end


% Write completion message
fprintf(fid, '\nAll models were enumerated successfully at %s.\n', datestr(now));

% Close the file
fclose(fid);
