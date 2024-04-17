function [model, In, Out] = nlsaModel_mlo(NLSA)
% Construct NLSA model for analysis of MLO CO2 data

experiment = experimentStr(NLSA);

switch experiment
    case 'CO2_197601-202211_emb48_l2_den'
        % NLSA parameters; in-sample data
        In.nN         = 0;          % nearest neighbors; defaults to max. value if 0
        In.nBQ        = 1;          % number of batches for query partition
        In.nBT        = 1;          % number of batches for test partition
        In.lDist      = 'l2';       % local distance
        In.tol        = 0;          % 0 distance threshold (for cone kernel)
        In.zeta       = 0.995;      % cone kernel parameter
        In.coneAlpha  = 0;          % velocity exponent in cone kernel
        In.nNS        = In.nN;      % nearest neighbors for symmetric distance
        In.diffOpType = 'gl_mb_bs'; % diffusion operator type
        In.epsilon    = 1;          % kernel bandwidth parameter
        In.epsilonB   = 2;          % kernel bandwidth base
        In.epsilonE   = [-40 40]; % kernel bandwidth exponents
        In.nEpsilon   = 200;        % number of exponents for bandwidth tuning
        In.alpha      = 0.5;        % diffusion maps normalization
        In.nPhi       = 201;         % diffusion eigenfunctions to compute
        In.nPhiPrj    = In.nPhi;    % eigenfunctions to project the data
        In.idxPhiRec  = [2 3];    % eigenfunctions for reconstruction
        In.idxPhiSVD  = 1 : 1;        % eigenfunctions for linear mapping
        In.idxVTRec   = 1 : 1;        % SVD termporal patterns for reconstruction

        % NLSA parameters, kernel density estimation
        In.denType     = 'vb';          % density estimation type
        In.denND       = 5;             % manifold dimension
        In.denLDist    = 'l2';          % local distance function
        In.denBeta     = -1 / In.denND; % density exponent
        In.denNN       = 50;            % nearest neighbors
        In.denZeta     = 0;             % cone kernel parameter
        In.denConeAlpha= 0;             % cone kernel velocity exponent
        In.denEpsilon  = 1;             % kernel bandwidth
        In.denEpsilonB = 2;             % kernel bandwidth base
        In.denEpsilonE = [-40 40];      % kernel bandwidth exponents
        In.denNEpsilon = 200;      % number of exponents for bandwidth tuning

        % Koopman generator parameters; in-sample data
        In.koopmanOpType = 'diff';     % Koopman generator type
        In.koopmanFDType  = 'central'; % finite-difference type
        In.koopmanFDOrder = 4;         % finite-difference order
        In.koopmanDt      = 1;         % sampling interval (in months)
        In.koopmanAntisym = true;      % enforce antisymmetrization
        In.koopmanEpsilon = 1.0E-4;      % regularization parameter
        In.koopmanRegType = 'inv';     % regularization type
        In.idxPhiKoopman  = 1 : 201;   % diffusion eigenfunctions used as basis
        In.nPhiKoopman    = numel(In.idxPhiKoopman); % eigenfunctions to compute
        In.nKoopmanPrj    = In.nPhiKoopman; % eigenfunctions to project the data
        In.idxKoopmanRec  = [ 2 3 ]; % eigenfunctions for recontruction

    otherwise
        error('Invalid experiment.')
end


%% PREPARE SOURCE AND TARGET COMPONENTS
% Use dummy lat/lon limits for station data

% Delay-embedding/finite-difference parameters; in-sample source data
In.Src(1).idxE      = 1 : NLSA.embWindow;     % delay-embedding indices
In.Src(1).nXB       = 0;          % samples before main interval
In.Src(1).nXA       = 0;          % samples after main interval
In.Src(1).fdOrder   = 0;          % finite-difference order
In.Src(1).fdType    = 'central';  % finite-difference type
In.Src(1).embFormat = 'overlap';  % storage format

% Delay-embedding/finite-difference parameters; in-sample target data
In.Trg(1).idxE      = 1 : 1;     % delay-embedding indices
In.Trg(1).nXB       = 0;          % samples before main interval
In.Trg(1).nXA       = 0;          % samples after main interval
In.Trg(1).fdOrder   = 0;          % finite-difference order
In.Trg(1).fdType    = 'central';  % finite-difference type
In.Trg(1).embFormat = 'overlap';  % storage format

if ~iscell(NLSA.var)
    NLSA.var = {NLSA.var};
end
nC = numel(NLSA.var);
for iC = nC : -1 : 1
    In.Src(iC).field = NLSA.var{iC};
    In.Src(iC).xLim = [1 1];
    In.Src(iC).yLim = [1 1];
    In.Trg(iC).field = NLSA.var{iC};
    In.Trg(iC).xLim = [1 1];
    In.Trg(iC).yLim = [1 1];
end


%% SERIAL DATE NUMBERS FOR IN-SAMPLE DATA
tNum1 = datenum(NLSA.tLim{1}, 'yyyymm');
tNum2 = datenum(NLSA.tLim{2}, 'yyyymm');
In.Res.experiment = 'mlo';
In.Res.tFormat = 'yyyymm';
In.Res.tLim = NLSA.tLim;
In.Res.nB = 1; % partition batches
In.Res.nBRec = 1; % batches for reconstructed data
nS = months(tNum1, tNum2) + 1;
In.Res.tNum = datemnth(tNum1, 0 : nS - 1);


%% SERIAL DATE NUMBERS FOR OUT-OF-SAMPLE DATA
ifOse = exist('Out', 'var');
if ifOse
    % TODO
end


%% CONSTRUCT NLSA MODEL
if ifOse
    args = {In Out};
else
    args = {In};
end
[model, In, Out] = climateNLSAModel(args{:});
