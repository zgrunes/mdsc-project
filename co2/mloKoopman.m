% KOOPMAN ANALYSIS OF MLO CO2 DATA
%
% Date range: 1973/08 to 2023/12.
% First non-missing datum: 1974/05
%
% To retrieve NLSA basis functions: phi = getDiffusionEigenfunctions(model);
% To retrieve Koopman eigenfunctions: z = getKoopmanEigenfunctions(model);
% To retrieve eigenperiodds: T = getKoopmanEigenperiods(model);
% To read raw data: x = getData(model.srcComponent);
%
% Modified 2024-04-06

%% SCRIPT EXECUTION OPTIONS
ifTrainingData          = true; % read training data in appropriate format for NLSA
ifNLSA                  = true; % perform NLSA
ifKoopman               = true; % compute Koopman eigenfunctions
ifKoopmanReconstruction = false; % perform Koopman reconstruction


%% DATA & NLSA PARAMETERS
NLSA.tLim     = {'197601' '202211'}; % analysis inverval
NLSA.var       = 'CO2'; % input variables
NLSA.embWindow = 48; % delay embedding window (months)
NLSA.kernel    = 'l2'; % kernel type
NLSA.den       = true; % set true to use variable-bandwidth kernel

dataFunc  = @importData_mlo; % data import function
modelFunc = @nlsaModel_mlo; % function to build NLSA model

experiment = experimentStr(NLSA);          % string identifier
disp(['NLSA EXPERIMENT: ' experiment])


%% BATCH PROCESSING
iProc = 1; % current process
nProc = 1; % number of processes


%% DATA RECONSTRUCTION
Rec.iC = 1; % component to reconstruct
Rec.iB = 1; % batch to reconstruct


%% READ DATA
if ifTrainingData
    disp(['Reading training data using function ' func2str(dataFunc) '...'])
    t = tic;
    dataFunc(NLSA);
    toc(t)
end


%% BUILD NLSA MODEL
disp('Building NLSA model...')
t = tic;
model = modelFunc(NLSA);
toc(t)


%% PERFORM NLSA
% Output from each step is saved on disk.
if ifNLSA
    disp('Takens delay embedding for source data...'); t = tic;
    computeDelayEmbedding(model)
    toc(t)

    % The following step is needed only if we are using velocity-dependent
    % kernels.
    if isa(model.embComponent, 'nlsaEmbeddedComponent_xi')
        disp('Phase space velocity (time tendency of data)...'); t = tic;
        computeVelocity(model)
        toc(t)
    end

    % The following step is only needed if query partition is employed for
    % source data.
    if ~isempty(model.embComponentQ)
        disp('Forming query partition for source data...'); t = tic;
        computeDelayEmbeddingQ(model)
        toc(t)
    end

    % The following step is only needed if test partition is employed for source
    % data.
    if ~isempty(model.embComponentT)
        disp('Forming test partition for source data...'); t = tic;
        computeDelayEmbeddingT(model)
        toc(t)
    end

    % The following steps are needed only if we are using variable-bandwidth
    % kernels.
    if isa(model, 'nlsaModel_den')
        fprintf('Pairwise distances for density data, %i/%i...\n', ...
                  iProc, nProc);
        t = tic;
        computeDenPairwiseDistances(model, iProc, nProc)
        toc(t)

        disp('Distance normalization for kernel density estimation...');
        t = tic;
        computeDenBandwidthNormalization(model);
        toc(t)

        disp('Kernel bandwidth tuning for density estimation...'); t = tic;
        computeDenKernelDoubleSum(model);
        toc(t)

        disp('Kernel density estimation...'); t = tic;
        computeDensity(model);
        toc(t)

        % The next step is only needed if a query partition was used for the
        % density data.
        if ~isempty(model.denEmbComponentQ)
            disp('Density splitting...'); t = tic;
            computeDensitySplitting(model);
            toc(t)
        end

        disp('Takens delay embedding for density data...'); t = tic;
        computeDensityDelayEmbedding(model);
        toc(t)
        % The following step is only needed if query partition is employed for
        % density data.
        if ~isempty(model.denEmbComponentQ)
            disp('Forming query partition for density data...'); t = tic;
            computeDensityDelayEmbeddingQ(model)
            toc(t)
        end

        % The following step is only needed if test partition is employed for
        % density data.
        if ~isempty(model.denEmbComponentT)
            disp('Forming test partition for density data...'); t = tic;
            computeDensityDelayEmbeddingT(model)
            toc(t)
        end
    end

    fprintf('Pairwise distances (%i/%i)...\n', iProc, nProc); t = tic;
    computePairwiseDistances(model, iProc, nProc)
    toc(t)

    disp('Distance symmetrization...'); t = tic;
    symmetrizeDistances(model)
    toc(t)

    disp('Kernel bandwidth tuning...'); t = tic;
    computeKernelDoubleSum(model)
    toc(t)

    disp('Kernel eigenfunctions...'); t = tic;
    computeDiffusionEigenfunctions(model)
    toc(t)
end


%% COMPUTE EIGENFUNCTIONS OF KOOPMAN GENERATOR
if ifKoopman
    disp('Koopman eigenfunctions...'); t = tic;
    computeKoopmanEigenfunctions(model, 'ifLeftEigenfunctions', true)
    toc(t)
end


%% PERFORM KOOPMAN RECONSTRUCTION
if ifKoopmanReconstruction

    disp('Takens delay embedding...'); t = tic;
    computeTrgDelayEmbedding(model, Rec.iC)
    toc(t)

    disp('Projection of target data onto Koopman eigenfunctions...'); t = tic;
    computeKoopmanProjection(model, Rec.iC)
    toc(t)

    disp('Reconstruction of target data from Koopman eigenfunctions...');
    t = tic;
    computeKoopmanReconstruction(model, Rec.iB, Rec.iC)
    toc(t)
end

%% Initialize the Koopman Analysis Script
% Make sure to compute the eigenfunctions and their periods
if ifKoopman
    disp('Retrieving Koopman eigenfunctions and eigenperiods...');
    z = getKoopmanEigenfunctions(model);  % Eigenfunctions
    T = getKoopmanEigenperiods(model);    % Eigenperiods
else
    error('Koopman analysis is not enabled. Set ifKoopman to true.');
end

%% Analyze Eigenvalues and Eigenperiods to Identify Trends and Oscillations

% Initialize parameters
nEigenfuncs = length(T);  % Number of eigenfunctions/eigenperiods
trendThreshold = 1e6;     % Threshold for what is considered a 'large' eigenperiod (can be adjusted)
oscillationThreshold = 12;  % Threshold for reasonable oscillation period in months (e.g., 12 months for annual)

% Preallocate arrays to store indices
trendIndices = [];
oscillationIndices = [];

% Loop through eigenperiods to categorize eigenfunctions
for i = 1:nEigenfuncs
    if isnan(T(i)) || T(i) > trendThreshold
        % Eigenfunctions with a very large or undefined period are trends
        trendIndices = [trendIndices, i];
    elseif T(i) > 0 && T(i) <= oscillationThreshold
        % Eigenfunctions with a defined, reasonable period are oscillations
        oscillationIndices = [oscillationIndices, i];
    end
end

%% Output Results
disp('Identified Trend-related Eigenfunctions:');
disp(trendIndices);

disp('Identified Oscillation-related Eigenfunctions:');
disp(oscillationIndices);

%% Optional: Visualization of Trend and Oscillation Eigenfunctions
figure;
subplot(1,2,1);
plot(real(z(:, trendIndices)), 'LineWidth', 2);
title('Trend Eigenfunctions');
xlabel('Time (months)');
ylabel('Eigenfunction Amplitude');

subplot(1,2,2);
plot(real(z(:, oscillationIndices)), 'LineWidth', 2);
title('Oscillation Eigenfunctions');
xlabel('Time (months)');
ylabel('Eigenfunction Amplitude');

sgtitle('Visualization of Identified Koopman Eigenfunctions');

