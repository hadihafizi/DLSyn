function determineSigSyn(irec,itime,triadChunk,dataToUse)
% This code will determine significant synergy values by computing synergy
% values from jittered spike trains and comparing those to observed synergy
% values. Jittering will be repeated 1000 times or until we have reached 10
% values which are larger than the observed value. This will allow us to
% enforce a p<0.01 significance for observed synergy values.

% There are many steps in this process which must be repeated for each
% triad, in each recording (25 total), at two timescales. These steps are
% listed and described below:
% Step 1: Bin the (triad) data: place spike times in bins according to the
% timescale of interest
% Step 2: Jitter the (transmitter) data: use a Monte Carlo approach which randomly
% selects values from a normal distribution to jitter spike times
% in a time window proportional to the timescale.
% Step 3: State the data: classify the data into states 1 or 2. State 1
% being not firing (time bins with no spikes) and state 2 being firing
% (time bins with spikes).
% Step 4: Compute synergy values using instinfo.m. For each jittered triad
% (trans only) compare to observed synergy value. Continue this process
% until we get through all 1000 jitterings or until we have generated 10
% jittered values that are larger than observed values. If this happens,
% stop the process and move on to the next triad (assign p-value > 0.01).

% Sam Faber
% Updated 02/24/2017

% directory for saving P-values
saveDir = '/N/u/samfaber/Karst/SynergyProject/results';

% point to spike time data
origDataDir = '/N/u/samfaber/Karst/Dropbox/MATLAB/Info Theory Network Project/TimmeCRCNSSubmissionVer1/Data';


if itime == 2
    binSize = 1.6;
else
    binSize = 3.5;
end

load(fullfile(origDataDir,['DataSet',num2str(irec),'.mat']));

% preallocate synery p value cell array
pValues = nan(size(dataToUse,1),1);

% Step 1: Bin all receiver data - use Nick's code parameter settings
maxLead = 0;
maxLag = data.recordinglength;
binEdges = (-floor(maxLead/binSize):floor(maxLag/binSize))*binSize; %data.recordinglength;
binnedRecData = zeros(data.nNeurons,(length(binEdges)-1));

receiverInds = unique(dataToUse(:,1));

for ireceiver = 1:length(receiverInds)
    dataBins = discretize(data.spikes{receiverInds(ireceiver)},binEdges,'IncludedEdge','Right'); % gives indices of time bins that contain spikes
    [nSpikes,perBin] = hist(dataBins,unique(dataBins)); % gives the number of spikes and the bins in which they occur
    binnedRecData(receiverInds(ireceiver),perBin)=nSpikes;
end
binnedRecData = uint8(binnedRecData);

% Step 2: Jitter transmitter data up to 1000 times
for itriad = 1:length(dataToUse); % for each triad
    
    % NOTE: instead of jittering transmitter spike trains 1000 times, we
    % will jitter each transmitter spike train 46 times and use
    % unique combinations of jittered transmitter time series to
    % compute synergy. We will use 1000 unique combinations which
    % would be equivalent to doing 1000 'jitterations.' This should
    % cut down on computation time.
    
    
    % get inds to triad transmitters
    transInds = dataToUse(itriad,2:3);
    nJitters = 46; 
    
    % first, jitter all transmitter data 46 times
    ASDF = cell(nJitters*length(transInds)+2,1);
    
    for ijitter = 1:(nJitters*length(transInds)) % n jitters X n transmitters
        
        if mod(ijitter,2)==0
            ASDF{ijitter} = data.spikes{transInds(2)}/binSize; %binnedData(transInds(itrans),:);
        else
            ASDF{ijitter} = data.spikes{transInds(1)}/binSize;
        end
    end
    ASDF{end-1} = binSize;
    ASDF{end} = [nJitters*length(transInds), data.recordinglength/binSize];
    
    
    % Set jittering parameters (according to Nick)
    Mode = 'Uniform';
    Value = 3;% number of bins to right or left that spike could be moved
    Preserve = 2; % preserves number of spikes in past states
    
    % perform jittering - this structure contains all the
    % jittered data that we need for computing synergy
    [ RandASDF ] = ASDFJitter(ASDF,Mode,Value,0.05,Preserve);
    
    % This would be a good place to perform checks to ensure that
    % jittering is taking place and that spikes are not being lost
    
    
    % Bin the jittered transmitter data
    maxLead = 0;
    maxLag = data.recordinglength;
    binEdges = (-floor(maxLead/binSize):floor(maxLag/binSize))*binSize; %data.recordinglength;
    binnedTransData = zeros(nJitters*length(transInds),(length(binEdges)-1));
    
    for itrans = 1:nJitters*length(transInds)
        dataBins = discretize(RandASDF{itrans},binEdges,'IncludedEdge','Right'); % gives indices of time bins that contain spikes
        [nSpikes,perBin] = hist(dataBins,unique(dataBins)); % gives the number of spikes and the bins in which they occur
        binnedTransData(itrans,perBin)=nSpikes;
    end
    
    % change array precision to cut down on memory usage
    binnedTransData = uint8(binnedTransData);
    transB = nan(nJitters,data.recordinglength/binSize);
    transA = nan(nJitters,data.recordinglength/binSize);
    
    % separate out combined binned transmitter data into two, separate
    % transmitter datasets
    for ii = 1:nJitters
        if mod(ii,2)==0
            transB(ii/2,:) = binnedTransData(ii,:);
        else
            transA((ii/2+(1/2)),:) = binnedTransData(ii,:);
        end
    end
    
    % we have all the jittered (separated) transmitter data
    transA = uint8(transA);
    transB = uint8(transB);
    
    % Now, we will prepare to begin synergy calculations
    
    % create unique combinations of transmitters
    combs = nchoosek(1:46,2);
    indsToJitteredTrans = combs(randperm(1000),:);
    
    % set the number of our first iteration for the while loop
    irepeat = 1;
   
    % when we start with a new triad, the number of big jitters is
    % zero
    nBigJitters = 0;
    
    
    % Step 3: State the jittered data (for the transmitters) and the
    % original data (for the receiver)
    
    %  isolate triad receiver data
    recData = binnedRecData(dataToUse(itriad,1),:);
    
    % preallocate dataraster
    DataRaster = nan(3,data.recordinglength/binSize);
    StatesRaster = nan(3,data.recordinglength/binSize);
    
    % compute synergy up to 1000 times
    while irepeat <= 1000 % number of jitter attempts
        
        % isolate triad transmitter data
        transAData = transA(indsToJitteredTrans(irepeat,1),:);
        transBData = transB(indsToJitteredTrans(irepeat,2),:);
        
        % get DataRaster
        DataRaster = [recData;transAData;transBData];
        StatesRaster = DataRaster;
        
        % state triad data
        DataRaster(DataRaster >= 1) = 2;
        DataRaster(DataRaster == 0) = 1;
        StatesRaster(:,1:(end - 1)) = StatesRaster(:,1:(end - 1)) + StatesRaster(:,2:end); 
        
        % Splice past time bins together
        StatesRaster(StatesRaster >= 1) = 2;
        StatesRaster(StatesRaster == 0) = 1;
        
        
        % Step 4: Compute synergy
        
        % compute synergy
        Method = '3TE';
        VariableIDs = {1, 1, 4 ;... % Receiving variable in the future (4 results from the delay bin and from the past states being 2 bins instead of 1)
            2, 1, 1 ;... % Receiving variable in the past
            2, 2, 1 ;  % Transmitting variable in the past
            2, 3, 1 }; % Other Transmitting variable in the past
        [InfoResult] = instinfo_SamVer({DataRaster;StatesRaster}, Method, VariableIDs);
        jitteredSynVal = InfoResult(4);
        
        % keep track of jitters
        if (jitteredSynVal - dataToUse(itriad,4)) > 100*eps % check if shuffled syn > observed syn value
            nBigJitters = nBigJitters + 1;
        end
        
        % if we surpass the significance threshold, stop
        if nBigJitters == 10
            pValues(itriad) = 0.01;
            irepeat = 1000 + 1; % don't compute any more synergy values for this triad, we are done
        else irepeat = irepeat + 1;
        end
        if irepeat == 1000 && nBigJitters < 10
            pValues(itriad) = nBigJitters/1000; % record the p-value
            irepeat = irepeat + 1; % we are done computing synergy for this triad
        end
    end
end

% save synergy p-values
save(fullfile(saveDir,['synergyPVals_rec_',num2str(irec),'_time_',num2str(itime),...
    '_triads_',num2str(triadChunk(1)),'_to_',num2str(triadChunk(2)) '.mat']),'pValues');




