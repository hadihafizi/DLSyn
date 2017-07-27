

clear all

dataset = {'139','151','152','163','165','168','174'};

for i = 1:length(dataset)
    asdfStruct = load(['asdf_rest',dataset{i},'_TimeCorOneMSGap.mat']);
%     load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat']);
%     load([dataset{i},'/wgts_1_16ms.mat'],'wgt');
%     Adj_mat = PDF_cor.*wgt; Adj_mat(Adj_mat~=0) = 1;
    
    asdf = ASDFChangeBinning(asdfStruct.asdf_raw,1); %clear asdf_raw
    asdf{end-1} = 1;
    
    spikes = asdf(1:end-2); recordingLength = asdf{end}(2);
    nNeurons = asdf{end}(1); binSize = 1.6; delayBin = 4;
%     triads = findPossibleTriads(Adj_mat);
    
    tic
    [TE] = computeTE(spikes, recordingLength, nNeurons, binSize, delayBin);
    toc
    varNames = {'TE'};
    varVals = {TE};
    parsave(['TE_',dataset{i},'.mat'],varNames,varVals);
    
    Adj_mat = TE; Adj_mat(Adj_mat ~= 0) = 1;
    triads = findPossibleTriads(Adj_mat);
    
    tic
    [Synergy, RecEnt] = computeSynergy(spikes, recordingLength, nNeurons, binSize, triads, delayBin);
    toc
    
    varNames = {'TE','triads','Synergy','RecEnt'};
    varVals = {TE,triads,Synergy,RecEnt};
    parsave(['TESynergyTriads_',dataset{i},'.mat'],varNames,varVals);
    
%     clearvars -except i dataset
    
end