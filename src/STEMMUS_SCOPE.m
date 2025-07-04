%% STEMMUS-SCOPE.m (script)
%     STEMMUS-SCOPE is a model for Integrated modeling of canopy photosynthesis, fluorescence,
%     and the transfer of energy, mass, and momentum in the soil-plant-atmosphere continuum
%     Copyright (C) 2021  Yunfei Wang, Lianyu Yu, Yijian Zeng, Christiaan Van der Tol, Bob Su
%     Contact: y.wang-3@utwente.nl; l.yu@utwente.nl; y.zeng@utwente.nl; c.vandertol@utwente.nl; z.su@utwente.nl
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public    License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Load in required Octave packages if STEMMUS-SCOPE is being run in Octave:
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    disp('Loading Octave packages...');
    pkg load statistics io;
end

% set CFG to a path if it is not defined
if exist('CFG', 'var') == 0
    CFG = "/home/khanalp/STEMMUS_SCOPE_model/STEMMUS_SCOPE_old/STEMMUS_SCOPE/ICOS_sites/IT-Lav/input/IT-Lav_2025-07-01-1331/IT-Lav_2025-07-01-1331_config.txt";
end

% set runMode to "full" if it is not defined
if exist('runMode', 'var') == 0
    runMode = "full";
    bmiMode = "none";
end

% Initialization routine
start_time = clock;
if strcmp(bmiMode, "initialize") || strcmp(runMode, "full")
    % Read the configPath file. Due to using MATLAB compiler, we cannot use run(CFG)
    disp (['Reading config from ', CFG]);
    [InputPath, OutputPath, InitialConditionPath, FullCSVfiles] = io.read_config(CFG);

    % Prepare forcing and soil data
    [SiteProperties, SoilProperties, TimeProperties] = io.prepareInputData(InputPath);

    % Load model settings: replacing "run Constants"
    ModelSettings = io.getModelSettings();

    % Load soil layers settings
    soilLayersFile = [InputPath, 'input_soilLayThick.csv'];
    if isfile(soilLayersFile)
        SoilLayerSettings = io.readSoilLayerSettings(soilLayersFile);
        % Override the default settings
        fields = fieldnames(SoilLayerSettings);
        for i = 1:numel(fields)
            ModelSettings.(fields{i}) = SoilLayerSettings.(fields{i});
        end
    end

    % Calculate other ModelSettings
    ModelSettings.NN = ModelSettings.NL + 1; % Number of nodes;
    ModelSettings.mN = ModelSettings.NN + 1;
    ModelSettings.mL = ModelSettings.NL + 1; % Number of elements. Preventing the exceeds of arrays size;
    ModelSettings.nD = 2;

    % Load groundwater settings
    GroundwaterSettings = groundwater.initializeGroundwaterSettings();
    GroundwaterSettings.soilThick = io.calculateSoilLayerThickness(ModelSettings);

    % load forcing data
    ForcingData = io.loadForcingData(InputPath, TimeProperties, SoilProperties, ModelSettings.Tot_Depth, GroundwaterSettings);

    % Get initial values
    InitialValues = init.defineInitialValues(TimeProperties.Dur_tot, ModelSettings);
    T_CRIT = InitialValues.T_CRIT;
    TT_CRIT = InitialValues.TT_CRIT;
    hOLD = InitialValues.hOLD;
    TOLD = InitialValues.TOLD;
    SAVE = InitialValues.SAVE;
    P_gOLD = InitialValues.P_gOLD;
    gwfluxes = struct();
    Evap = [];
    Trap = [];
    RWUs = [];
    RWUg = [];

    %% 1. define Constants
    Constants = io.define_constants();

    RTB = 1000; % initial root total biomass (g m-2)
    % Rl used in ebal
    [Rl, Ztot] = Initial_root_biomass(RTB, ModelSettings.DeltZ_R, ModelSettings.rroot, ModelSettings.ML, SiteProperties.landcoverClass(1));

    %% 2. simulation options
    path_input = InputPath;  % path of all inputs
    path_of_code = cd;

    parameter_file = {'input_data.xlsx'};
    options = io.readStructFromExcel([path_input char(parameter_file)], 'options', 2, 1);

    if options.simulation > 2 || options.simulation < 0
        fprintf('\n simulation option should be between 0 and 2 \r');
        return
    end

    %% 3. file names
    [dummy, X] = xlsread([path_input char(parameter_file)], 'filenames');
    j = find(~strcmp(X(:, 2), {''}));
    X = X(j, 1:end);

    F = struct('FileID', {'Simulation_Name', 'soil_file', 'leaf_file', 'atmos_file'...
                          'Dataset_dir', 't_file', 'year_file', 'Rin_file', 'Rli_file' ...
                          'p_file', 'Ta_file', 'ea_file', 'u_file', 'CO2_file', 'z_file', 'tts_file' ...
                          'LAI_file', 'hc_file', 'SMC_file', 'Vcmax_file', 'Cab_file', 'LIDF_file'});
    for i = 1:length(F)
        k = find(strcmp(F(i).FileID, strtok(X(:, 1))));
        if ~isempty(k)
            F(i).FileName = strtok(X(k, 2));
        end
    end

    %% 4. input data
    [N, X] = xlsread([path_input char(parameter_file)], 'inputdata', '');
    X  = X(9:end, 1);

    % Create a structure holding Scope parameters
    useXLSX = 1; % set it to 1 or 0, the current stemmus-scope does not support 0
    [ScopeParameters, options] = parameters.loadParameters(options, useXLSX, X, F, N);

    % Define the location information
    ScopeParameters.LAT = SiteProperties.latitude; % latitude
    ScopeParameters.LON = SiteProperties.longitude; % longitude
    ScopeParameters.BSMlat = SiteProperties.latitude; % latitude of BSM model
    ScopeParameters.BSMlon = SiteProperties.longitude; % longitude of BSM model
    ScopeParameters.z =  SiteProperties.reference_height;   % reference height
    ScopeParameters.hc =  SiteProperties.canopy_height;  % canopy height
    ScopeParameters.Tyear = mean(ForcingData.Ta_msr); % calculate mean air temperature

    % calculate the time zone based on longitude
    ScopeParameters.timezn = helpers.calculateTimeZone(SiteProperties.longitude);

    % Input T parameters for different vegetation type
    [ScopeParameters] = parameters.setTempParameters(ScopeParameters, SiteProperties.sitename, SiteProperties.landcoverClass);

    %% 5. Declare paths
    path_input      = InputPath;          % path of all inputs
    path_output     = OutputPath;
    %% 6. Numerical parameters (iteration stops etc)
    iter.maxit           = 100;                          %                   maximum number of iterations
    iter.maxEBer         = 5;                            % [W m-2]            maximum accepted error in energy bal.
    iter.Wc              = 1;                         %                   Weight coefficient for iterative calculation of Tc

    %% 7. Load spectral data for leaf and soil
    load([path_input, 'fluspect_parameters/', char(F(3).FileName)]);
    rsfile      = load([path_input, 'soil_spectrum/', char(F(2).FileName)]);        % file with soil reflectance spectra

    %% 8. Load directional data from a file
    directional = struct;
    if options.calc_directional
        anglesfile          = load([path_input, 'directional/brdf_angles2.dat']); %     Multiple observation angles in case of BRDF calculation
        directional.tto     = anglesfile(:, 1);              % [deg]             Observation zenith Angles for calcbrdf
        directional.psi     = anglesfile(:, 2);              % [deg]             Observation zenith Angles for calcbrdf
        directional.noa     = length(directional.tto);      %                   Number of Observation Angles
    end

    %% 9. Define canopy structure
    canopy.nlayers  = 30;
    nl              = canopy.nlayers;
    canopy.x        = (-1 / nl:-1 / nl:-1)';         % a column vector
    canopy.xl       = [0; canopy.x];                 % add top level
    canopy.nlincl   = 13;
    canopy.nlazi    = 36;
    canopy.litab    = [5:10:75 81:2:89]';   % a column, never change the angles unless 'ladgen' is also adapted
    canopy.lazitab  = (5:10:355);           % a row

    %% 10. Define spectral regions
    [spectral] = io.define_bands();

    wlS  = spectral.wlS;    % SCOPE 1.40 definition
    wlP  = spectral.wlP;    % PROSPECT (fluspect) range
    wlT  = spectral.wlT;    % Thermal range
    wlF  = spectral.wlF;    % Fluorescence range

    I01  = find(wlS < min(wlF));   % zero-fill ranges for fluorescence
    I02  = find(wlS > max(wlF));
    N01  = length(I01);
    N02  = length(I02);

    nwlP = length(wlP);
    nwlT = length(wlT);

    nwlS = length(wlS);

    spectral.IwlP = 1:nwlP;
    spectral.IwlT = nwlP + 1:nwlP + nwlT;
    spectral.IwlF = (640:850) - 399;

    [rho, tau, rs] = deal(zeros(nwlP + nwlT, 1));

    %% 11. load time series data
    ScopeParametersNames = fieldnames(ScopeParameters);
    if options.simulation == 1
        vi = ones(length(ScopeParametersNames), 1);
        [soil, leafbio, canopy, meteo, angles, xyt]  = io.select_input(ScopeParameters, InitialValues.Theta_LL, vi, canopy, options, SiteProperties, SoilProperties);
        [ScopeParameters, xyt, canopy]  = io.loadTimeSeries(ScopeParameters, leafbio, soil, canopy, meteo, F, xyt, path_input, options, SiteProperties.landcoverClass);
    else
        soil = struct;
    end

    %% 12. preparations
    if options.simulation == 1
        diff_tmin           =   abs(xyt.t - xyt.startDOY);
        diff_tmax           =   abs(xyt.t - xyt.endDOY);
        I_tmin              =   find(min(diff_tmin) == diff_tmin);
        I_tmax              =   find(min(diff_tmax) == diff_tmax);
        if options.soil_heat_method < 2
            meteo.Ta = ForcingData.Ta_msr(1);
            soil.Tsold = meteo.Ta * ones(12, 2);
        end
    end
    ScopeParametersNames = fieldnames(ScopeParameters);
    nvars = length(ScopeParametersNames);
    vmax = ones(nvars, 1);
    for i = 1:nvars
        name = ScopeParametersNames{i};
        vmax(i) = length(ScopeParameters.(name));
    end
    vmax(27, 1) = 1; % these are Tparam and LIDFb
    vi      = ones(nvars, 1);
    switch options.simulation
        case 0
            telmax = max(vmax);
            [xyt.t, xyt.year] = deal(zeros(telmax, 1));
        case 1
            telmax = size(xyt.t, 1);
        case 2
            telmax  = prod(double(vmax));
            [xyt.t, xyt.year] = deal(zeros(telmax, 1));
    end
    [rad, thermal, fluxes] = io.initialize_output_structures(spectral);
    atmfile     = [path_input 'radiationdata/' char(F(4).FileName(1))];
    atmo.M      = helpers.aggreg(atmfile, spectral.SCOPEspec);

    %% 13. create output files and
    [Output_dir, fnames] = io.create_output_files_binary(parameter_file, SiteProperties.sitename, path_of_code, path_input, path_output, spectral, options);

    %% Initialize Temperature, Matric potential and soil air pressure.
    % Define soil variables for StartInit
    VanGenuchten = init.setVanGenuchtenParameters(SoilProperties);
    SoilVariables = init.defineSoilVariables(InitialValues, SoilProperties, VanGenuchten);

    % Add initial soil moisture and soil temperature
    [SoilInitialValues, BtmX, BtmT, Tss] = io.loadSoilInitialValues(InputPath, TimeProperties, SoilProperties, ForcingData, ModelSettings);
    SoilVariables.InitialValues = SoilInitialValues;
    SoilVariables.BtmX = BtmX;
    SoilVariables.BtmT = BtmT;
    SoilVariables.Tss = Tss;

    % Run StartInit
    [SoilVariables, VanGenuchten, ThermalConductivity] = StartInit(SoilVariables, SoilProperties, VanGenuchten, ModelSettings);

    % Set SoilVariables that are used in the loop
    T = SoilVariables.T;
    h = SoilVariables.h;
    TT = SoilVariables.TT;
    h_frez = SoilVariables.h_frez;
    hh = SoilVariables.hh;
    hh_frez = SoilVariables.hh_frez;

    % get soil constants
    SoilConstants = io.getSoilConstants();

    %% The boundary condition information settings
    BoundaryCondition = init.setBoundaryCondition(SoilVariables, ForcingData, SiteProperties.landcoverClass(1), ModelSettings);
    DSTOR = BoundaryCondition.DSTOR;
    DSTOR0 = BoundaryCondition.DSTOR0;
    RS = BoundaryCondition.RS;
    DSTMAX = BoundaryCondition.DSTMAX;
    IRPT1 = BoundaryCondition.IRPT1;
    IRPT2 = BoundaryCondition.IRPT2;

    %% 14. Initialize matrices
    calculate = 1;
    TIMEOLD = 0;
    SAVEhh_frez = zeros(ModelSettings.NN + 1, 1);
    FCHK = zeros(1, ModelSettings.NN);
    KCHK = zeros(1, ModelSettings.NN);
    hCHK = zeros(1, ModelSettings.NN);
    TIMELAST = 0;

    % the start of simulation period is from 0mins, while the input data start from 30mins.
    TIME = 0;  % Time of simulation released;
    Delt_t = TimeProperties.DELT;  % Duration of time step [Unit of second]
    KT = ModelSettings.KT;
    kk = 0;
    TimeStep = [];
    TEND = TIME + TimeProperties.DELT * TimeProperties.Dur_tot; % Time to be reached at the end of simulation period
    Delt_t0 = Delt_t; % Duration of last time step
    TOLD_CRIT = [];

    % for soil moisture and temperature outputs
    monitorDepthTemperature = ModelSettings.NL:-1:1;
    monitorDepthSoilMoisture = ModelSettings.NL:-1:1;
    Sim_Theta_U = [];
    Sim_Temp = [];

    % for matric potential and fluxes outputs, added by Mostafa
    Sim_hh = [];
    Sim_qlh = [];
    Sim_qlt = [];
    Sim_qvh = [];
    Sim_qvt = [];
    Sim_qla = [];
    Sim_qva = [];

    % Srt, root water uptake;
    Srt = InitialValues.Srt;  % will be updated!
    P_gg = InitialValues.P_gg;  % will be updated!
    P_g = InitialValues.P_g;  % will be updated!

    %% soil layer information
    SoilLayer.thickness = ModelSettings.DeltZ_R;
    SoilLayer.depth = Ztot';  % used in Initial_root_biomass

    % NOTE: workspace will be saved, this code block is after the update step.
    % (this is to not repeat the save-workspace code).
    disp('Finished model initialization');
end

% If the runMode is update, retrieve the (possibly) updated state.
% The state can be modified by the STEMMUS_SCOPE BMI. See PyStemmusScope.
if strcmp(bmiMode, 'update')
    load([OutputPath, 'STEMMUS_SCOPE_state.mat']); % Load the workspace to be able to (continue) running the model

    if KT + 1 >= TimeProperties.Dur_tot
        bmiMode = 'none';  % Ensure the model does not try to update.
        disp("Finished running the model. Updating won't do anything!");
    else
        endTime = KT + 1;
    end
elseif strcmp(runMode, 'full')
    endTime = TimeProperties.Dur_tot;
    disp('The calculations start now');
end

% Actually run the model
if strcmp(bmiMode, 'update') || strcmp(runMode, 'full')

    if GroundwaterSettings.GroundwaterCoupling == 1  % Groundwater coupling is enabled
        BoundaryCondition.NBChB = 1;

        % Update GroundwaterSettings.headBotmLayer and GroundwaterSettings.tempBotm, from MODFLOW through BMI
        GroundwaterSettings.gw_Dep = groundwater.calculateGroundWaterDepth(GroundwaterSettings.topLevel, GroundwaterSettings.headBotmLayer, ModelSettings.Tot_Depth);

        % Update Dunnian runoff and ForcingData.Precip_msr
        [ForcingData.runoffDunnian, ForcingData.Precip_msr] = groundwater.updateDunnianRunoff(ForcingData.Precip_msr, GroundwaterSettings.gw_Dep);

        % Calculate the index of the bottom layer level
        [GroundwaterSettings.indxBotmLayer, GroundwaterSettings.indxBotmLayer_R] = groundwater.calculateIndexBottomLayer(GroundwaterSettings.soilThick, GroundwaterSettings.gw_Dep, ModelSettings);

        % Check the position of the groundwater table
        if any(GroundwaterSettings.gw_Dep > ModelSettings.Tot_Depth) || any(isempty(GroundwaterSettings.indxBotmLayer))
            warning('Groundwater table is below the end of the soil column. Please enlarge the total soil thickness!');
            return
        end

        % Check if no available groundwater temperature data
        if isnan(GroundwaterSettings.tempBotm)
            GroundwaterSettings.tempBotm = SoilVariables.TT(GroundwaterSettings.indxBotmLayer);
        end
    end

    WaterPotential = struct('leaf', []);
    
    % Will do one timestep in "update mode", and run until the end if in "full run" mode.
    while KT < endTime
        
        KT = KT + 1;  % Counting Number of timesteps
        warning('Timestep = %d', KT)
        if KT > 1 && Delt_t > (TEND - TIME)
            Delt_t = TEND - TIME;  % If Delt_t is changed due to excessive change of state variables, the judgement of the last time step is excuted.
        end
        TIME = TIME + Delt_t;  % The time elapsed since start of simulation
        TimeStep(KT, 1) = Delt_t;
        SUMTIME(KT, 1) = TIME;
        Processing = TIME / TEND;
        NoTime(KT) = fix(SUMTIME(KT) / TimeProperties.DELT);
        if NoTime(KT) == 0
            k = 1;
        else
            k = NoTime(KT);
        end
        %%%%% Updating the state variables. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L_f = 0;  % ignore Freeze/Thaw, see issue 139
        TT_CRIT(ModelSettings.NN) = ModelSettings.T0; % unit K
        hOLD_frez = [];
        if IRPT1 == 0 && IRPT2 == 0 && SoilVariables.ISFT == 0
            for MN = 1:ModelSettings.NN
                hOLD_frez(MN) = h_frez(MN);
                h_frez(MN) = hh_frez(MN);
                TOLD_CRIT(MN) = T_CRIT(MN);
                T_CRIT(MN) = TT_CRIT(MN);

                hOLD(MN) = h(MN);
                h(MN) = hh(MN);
                DDhDZ(MN, KT) = InitialValues.DhDZ(MN);
                if ModelSettings.Thmrlefc == 1
                    TOLD(MN) = T(MN);
                    T(MN) = TT(MN);
                    TTT(MN, KT) = TT(MN);
                end
                if ModelSettings.Soilairefc == 1
                    P_gOLD(MN) = P_g(MN);
                    P_g(MN) = P_gg(MN);
                end
                if ModelSettings.rwuef == 1
                    SRT(MN, KT) = Srt(MN, 1);
                end
            end
            if options.simulation == 1
                vi(vmax > 1) = k;
            end
            if options.simulation == 0
                vi(vmax == telmax) = k;
            end
            [soil, leafbio, canopy, meteo, angles, xyt] = io.select_input(ScopeParameters, SoilVariables.Theta_LL, vi, canopy, options, xyt, soil);
            if options.simulation ~= 1
                fprintf('simulation %i ', k);
                fprintf('of %i \n', telmax);
            else
                calculate = 1;
            end

            if calculate
                iter.counter = 0;
                LIDF_file            = char(F(22).FileName);
                if  ~isempty(LIDF_file)
                    canopy.lidf     = dlmread([path_input, 'leafangles/', LIDF_file], '', 3, 0);
                else
                    canopy.lidf     = equations.leafangles(canopy.LIDFa, canopy.LIDFb);    % This is 'ladgen' in the original SAIL model,
                end

                leafbio.emis        = 1 - leafbio.rho_thermal - leafbio.tau_thermal;

                if options.calc_PSI
                    fversion = @fluspect_B_CX;
                else
                    fversion = @fluspect_B_CX_PSI_PSII_combined;
                end
                leafbio.V2Z = 0;
                leafopt  = fversion(spectral, leafbio, optipar);
                leafbio.V2Z = 1;
                leafoptZ = fversion(spectral, leafbio, optipar);
                IwlP     = spectral.IwlP;
                IwlT     = spectral.IwlT;
                rho(IwlP)  = leafopt.refl;
                tau(IwlP)  = leafopt.tran;
                rlast    = rho(nwlP);
                tlast    = tau(nwlP);

                if options.soilspectrum == 0
                    rs(IwlP) = rsfile(:, soil.spectrum + 1);
                else
                    soilemp.SMC   = 25;        % empirical parameter (fixed)
                    soilemp.film  = 0.015;     % empirical parameter (fixed)
                    rs(IwlP) = BSM(soil, optipar, soilemp);
                end
                rslast   = rs(nwlP);

                switch options.rt_thermal
                    case 0
                        rho(IwlT) = ones(nwlT, 1) * leafbio.rho_thermal;
                        tau(IwlT) = ones(nwlT, 1) * leafbio.tau_thermal;
                        rs(IwlT)  = ones(nwlT, 1) * soil.rs_thermal;
                    case 1
                        rho(IwlT) = ones(nwlT, 1) * rlast;
                        tau(IwlT) = ones(nwlT, 1) * tlast;
                        rs(IwlT)  = ones(nwlT, 1) * rslast;
                end
                leafopt.refl = rho;     % extended wavelength ranges are stored in structures
                leafopt.tran = tau;
                reflZ = leafopt.refl;
                tranZ = leafopt.tran;
                reflZ(1:300) = leafoptZ.refl(1:300);
                tranZ(1:300) = leafoptZ.tran(1:300);
                leafopt.reflZ = reflZ;
                leafopt.tranZ = tranZ;
                soil.refl    = rs;
                soil.Ts     = meteo.Ta * ones(2, 1);       % initial soil surface temperature

                if length(F(4).FileName) > 1 && options.simulation == 0
                    atmfile     = [path_input 'radiationdata/' char(F(4).FileName(k))];
                    atmo.M      = helpers.aggreg(atmfile, spectral.SCOPEspec);
                end
                atmo.Ta     = meteo.Ta;

                [rad, gap, profiles]   = RTMo(spectral, atmo, soil, leafopt, canopy, angles, meteo, rad, options);

                switch options.calc_ebal
                    case 1
                        [iter, fluxes, rad, thermal, profiles, soil, RWU, frac, WaterStressFactor, WaterPotential] ...
                            = ebal(iter, options, spectral, rad, gap,  ...
                                   leafopt, angles, meteo, soil, canopy, leafbio, xyt, k, profiles, Delt_t, ...
                                   Rl, SoilVariables, VanGenuchten, InitialValues, ModelSettings, GroundwaterSettings,WaterPotential,KT);
                        if options.calc_fluor
                            if options.calc_vert_profiles
                                [rad, profiles] = RTMf(spectral, rad, soil, leafopt, canopy, gap, angles, profiles);
                            else
                                [rad] = RTMf(spectral, rad, soil, leafopt, canopy, gap, angles, profiles);
                            end
                        end
                        if options.calc_xanthophyllabs
                            [rad] = RTMz(spectral, rad, soil, leafopt, canopy, gap, angles, profiles);
                        end

                        if options.calc_planck
                            rad         = RTMt_planck(spectral, rad, soil, leafopt, canopy, gap, angles, thermal.Tcu, thermal.Tch, thermal.Ts(2), thermal.Ts(1), 1);
                        end

                        if options.calc_directional
                            directional = calc_brdf(options, directional, spectral, angles, rad, atmo, soil, leafopt, canopy, meteo, profiles, thermal);
                        end

                    otherwise
                        Fc              = (1 - gap.Ps(1:end - 1))' / nl;      %           Matrix containing values for Ps of canopy
                        fluxes.aPAR     = canopy.LAI * (Fc * rad.Pnh        + equations.meanleaf(canopy, rad.Pnu, 'angles_and_layers', gap.Ps)); % net PAR leaves
                        fluxes.aPAR_Cab = canopy.LAI * (Fc * rad.Pnh_Cab    + equations.meanleaf(canopy, rad.Pnu_Cab, 'angles_and_layers', gap.Ps)); % net PAR leaves
                        [fluxes.aPAR_Wm2, fluxes.aPAR_Cab_eta] = deal(canopy.LAI * (Fc * rad.Rnh_PAR    + equations.meanleaf(canopy, rad.Rnu_PAR, 'angles_and_layers', gap.Ps))); % net PAR leaves
                        if options.calc_fluor
                            profiles.etah = ones(nl, 1);
                            profiles.etau = ones(13, 36, nl);
                            if options.calc_vert_profiles
                                [rad, profiles] = RTMf(spectral, rad, soil, leafopt, canopy, gap, angles, profiles);
                            else
                                [rad] = RTMf(spectral, rad, soil, leafopt, canopy, gap, angles, profiles);
                            end
                        end
                end

                if options.calc_fluor % total emitted fluorescence irradiance (excluding leaf and canopy re-absorption and scattering)
                    if options.calc_PSI
                        rad.Femtot = 1E3 * (leafbio.fqe(2) * optipar.phiII(spectral.IwlF) * fluxes.aPAR_Cab_eta + leafbio.fqe(1) * optipar.phiI(spectral.IwlF)  * fluxes.aPAR_Cab);
                    else
                        rad.Femtot = 1E3 * leafbio.fqe * optipar.phi(spectral.IwlF) * fluxes.aPAR_Cab_eta;
                    end
                end
            end
            if options.simulation == 2 && telmax > 1
                vi  = helpers.count(nvars, vi, vmax, 1);
            end
            if KT == 1
                if isreal(fluxes.Actot) && isreal(thermal.Tsave) && isreal(fluxes.lEstot) && isreal(fluxes.lEctot)
                    Acc = fluxes.Actot;
                    Tss = thermal.Tsave;
                else
                    Acc = 0;
                    fluxes.lEstot = 0;
                    fluxes.lEctot = 0;
                    Tss = ForcingData.Ta_msr(KT);
                end
            elseif NoTime(KT) > NoTime(KT - 1)
                if isreal(fluxes.Actot) && isreal(thermal.Tsave) && isreal(fluxes.lEstot) && isreal(fluxes.lEctot)
                    Acc = fluxes.Actot;
                    Tss = thermal.Tsave;
                else
                    Acc = 0;
                    fluxes.lEstot = 0;
                    fluxes.lEctot = 0;
                    Tss = ForcingData.Ta_msr(KT);
                end
            end

            DSTOR0 = DSTOR;
            BoundaryCondition.DSTOR0 = DSTOR0;

            if KT > 1
                SoilVariables.XWRE = updateWettingHistory(SoilVariables, VanGenuchten, ModelSettings);
            end

            for i = 1:ModelSettings.NL
                DDhDZ(i, KT) = InitialValues.DhDZ(i);
            end
        end
        if Delt_t ~= Delt_t0
            for MN = 1:ModelSettings.NN
                hh(MN) = h(MN) + (h(MN) - hOLD(MN)) * Delt_t / Delt_t0;
                TT(MN) = T(MN) + (T(MN) - TOLD(MN)) * Delt_t / Delt_t0;
            end
        end
        hSAVE = hh(ModelSettings.NN);
        TSAVE = TT(ModelSettings.NN);

        % set "hN" empty when the "if statement" below does not happen, see issue 98,
        % item 5
        hN = [];
        if BoundaryCondition.NBCh == 1
            hN = BoundaryCondition.BCh;
            hh(ModelSettings.NN) = hN;
            hSAVE = hN;
        elseif BoundaryCondition.NBCh == 2
            if BoundaryCondition.NBChh ~= 2
                if BoundaryCondition.BCh < 0
                    hN = DSTOR0;
                    hh(ModelSettings.NN) = hN;
                    hSAVE = hN;
                else
                    hN = -1e6;
                    hh(ModelSettings.NN) = hN;
                    hSAVE = hN;
                end
            end
        else
            if BoundaryCondition.NBChh ~= 2
                hN = DSTOR0;
                hh(ModelSettings.NN) = hN;
                hSAVE = hN;
            end
        end

        % Start the iteration procedure in a time step.
        SoilVariables.h = h;
        SoilVariables.hh = hh;
        SoilVariables.TT = TT;
        SoilVariables.T = T;
        SoilVariables.h_frez = h_frez;
        SoilVariables.Tss(KT) = Tss;

        for KIT = 1:ModelSettings.NIT
            [TT_CRIT, hh_frez] = HT_frez(SoilVariables.hh, ModelSettings.T0, Constants.g, L_f, SoilVariables.TT, ModelSettings.NN, SoilConstants.hd, ForcingData.Tmin);
            % update inputs for UpdateSoilWaterContent
            SoilVariables.TT_CRIT = TT_CRIT;
            SoilVariables.hh_frez = hh_frez;

            SoilVariables = UpdateSoilWaterContent(KIT, L_f, SoilVariables, VanGenuchten, ModelSettings);

            % Reset KL_T here. CondL_T script is replaced by this line
            % see issue 181, item 4
            KL_T = InitialValues.KL_T;

            [RHOV, DRHOVh, DRHOVT] = Density_V(SoilVariables.TT, SoilVariables.hh, Constants.g, Constants.Rv, ModelSettings.NN);

            TransportCoefficient = conductivity.calculateTransportCoefficient(InitialValues, SoilVariables, VanGenuchten, Delt_t, ModelSettings);
            W = TransportCoefficient.W;
            D_Ta = TransportCoefficient.D_Ta;

            [L] = Latent(SoilVariables.TT, ModelSettings.NN);
            % DRHODAt unused!
            [Xaa, XaT, Xah, DRHODAt, DRHODAz, RHODA] = Density_DA(SoilVariables.T, Constants.RDA, P_g, Constants.Rv, ModelSettings.DeltZ, SoilVariables.h, SoilVariables.hh, SoilVariables.TT, P_gg, Delt_t, ModelSettings.NL, ModelSettings.NN, DRHOVT, DRHOVh, RHOV);

            ThermalConductivityCapacity = conductivity.calculateThermalConductivityCapacity(InitialValues, ThermalConductivity, SoilVariables, VanGenuchten, ModelSettings, DRHOVT, L, RHOV);

            k_g = conductivity.calculateGasConductivity(InitialValues, TransportCoefficient, VanGenuchten, SoilVariables, ModelSettings);

            VaporVariables = conductivity.calculateVaporVariables(InitialValues, SoilVariables, VanGenuchten, ModelSettings, ThermalConductivityCapacity, SoilVariables.TT);

            GasDispersivity = conductivity.calculateGasDispersivity(InitialValues, SoilVariables, P_gg, k_g, ModelSettings);

            % Srt is both input and output
            [SoilVariables, HeatMatrices, HeatVariables, HBoundaryFlux, Rn_SOIL, Evap, Trap, r_a_SOIL, Srt, CHK, AVAIL0, RWUs, RWUg, ForcingData] = ...
                                  soilmoisture.solveSoilMoistureBalance(SoilVariables, InitialValues, ForcingData, VaporVariables, GasDispersivity, ...
                                                                        TimeProperties, SoilProperties, BoundaryCondition, Delt_t, RHOV, DRHOVh, ...
                                                                        DRHOVT, D_Ta, hN, RWU, fluxes, KT, KIT, hOLD, Srt, P_gg, ModelSettings, GroundwaterSettings);

            if BoundaryCondition.NBCh == 1
                DSTOR = 0;
                RS = 0;
            elseif BoundaryCondition.NBCh == 2
                AVAIL = -BoundaryCondition.BCh;
                EXCESS = (AVAIL + HBoundaryFlux.QMT) * Delt_t;
                if abs(EXCESS / Delt_t) <= 1e-10
                    EXCESS = 0;
                end
                DSTOR = min(EXCESS, DSTMAX);
                RS = (EXCESS - DSTOR) / Delt_t;
            else
                AVAIL = AVAIL0 - Evap;
                EXCESS = (AVAIL + HBoundaryFlux.QMT) * Delt_t; % (unit converstion from cm/sec to cm/30mins)
                if abs(EXCESS / Delt_t) <= 1e-10
                    EXCESS = 0;
                end
                DSTOR = min(EXCESS, DSTMAX); % Depth of depression storage at end of current time step
                % Next line is commented and Surface runoff is re-calculated using different approach (the following 3 lines)
                % RS(KT) = (EXCESS - DSTOR) / Delt_t; % surface runoff, (unit conversion from cm/30mins to cm/sec)
                ForcingData.runoff = ForcingData.runoffHort + ForcingData.runoffDunn; % total surface runoff (cm/sec)
            end

            if ModelSettings.Soilairefc == 1
                [AirVariabes, RHS, SAVE, P_gg] = dryair.solveDryAirEquations(SoilVariables, GasDispersivity, TransportCoefficient, InitialValues, VaporVariables, ...
                                                                             BoundaryCondition, ForcingData, P_gg, P_g, Xah, XaT, Xaa, RHODA, KT, Delt_t, ModelSettings, GroundwaterSettings);
            else
                AirVariabes.KLhBAR = InitialValues.KLhBAR;
                AirVariabes.KLTBAR = InitialValues.KLTBAR;
                AirVariabes.DDhDZ = DDhDZ;  % DDhDZ is not defined as InitialValues, see issue 100, item 3
                AirVariabes.DhDZ = InitialValues.DhDZ;
                AirVariabes.DTDZ = InitialValues.DTDZ;
                AirVariabes.Kaa = InitialValues.Kaa;
                AirVariabes.Vaa = InitialValues.Vaa;
                AirVariabes.QL = InitialValues.QL;
            end

            if ModelSettings.Thmrlefc == 1
                % CHK will be updated
                [RHS, SAVE, CHK, SoilVariables, EnergyVariables, gwfluxes] = energy.solveEnergyBalanceEquations(InitialValues, SoilVariables, HeatVariables, TransportCoefficient, ...
                                                                                                                AirVariabes, VaporVariables, GasDispersivity, ThermalConductivityCapacity, ...
                                                                                                                HBoundaryFlux, BoundaryCondition, ForcingData, DRHOVh, DRHOVT, KL_T, ...
                                                                                                                Xah, XaT, Xaa, Srt, L_f, RHOV, RHODA, DRHODAz, L, Delt_t, P_g, P_gg, ...
                                                                                                                TOLD, Evap, r_a_SOIL, Rn_SOIL, KT, CHK, ModelSettings, GroundwaterSettings);
            end

            if max(CHK) < 0.1
                break
            end
            hSAVE = SoilVariables.hh(ModelSettings.NN);
            TSAVE = SoilVariables.TT(ModelSettings.NN);
        end

        TIMEOLD = KT;
        KIT;
        KIT = 0;

        [TT_CRIT, hh_frez] = HT_frez(SoilVariables.hh, ModelSettings.T0, Constants.g, L_f, SoilVariables.TT, ModelSettings.NN, SoilConstants.hd, ForcingData.Tmin);
        % updates inputs for UpdateSoilWaterContent
        SoilVariables.TT_CRIT = TT_CRIT;
        SoilVariables.hh_frez = hh_frez;
        SoilVariables = UpdateSoilWaterContent(KIT, L_f, SoilVariables, VanGenuchten, ModelSettings);

        if IRPT1 == 0 && IRPT2 == 0
            if KT        % In case last time step is not convergent and needs to be repeated.
                for i = 1:ModelSettings.NL
                    hhh(i, KT) = SoilVariables.hh(i);
                    qlh(i, KT) = transpose(EnergyVariables.QL_h(i));
                    qlt(i, KT) = transpose(EnergyVariables.QL_T(i));
                    qla(i, KT) = transpose(EnergyVariables.QL_a(i));
                    qvh(i, KT) = transpose(EnergyVariables.QVH(i));
                    qvt(i, KT) = transpose(EnergyVariables.QVT(i));
                    qva(i, KT) = transpose(EnergyVariables.QVa(i));
                    qtot(i, KT) = transpose(EnergyVariables.Qtot(i));
                    for j = 1:ModelSettings.nD
                        Theta_LLL(i, j, KT) = SoilVariables.Theta_LL(i, j);
                        SoilVariables.Theta_L(i, j) = SoilVariables.Theta_LL(i, j);
                        Theta_UUU(i, j, KT) = SoilVariables.Theta_UU(i, j);
                        SoilVariables.Theta_U(i, j) = SoilVariables.Theta_UU(i, j);
                        Theta_III(i, j, KT) = SoilVariables.Theta_II(i, j);
                        SoilVariables.Theta_I(i, j) = SoilVariables.Theta_II(i, j);
                    end
                end

                % replace run ObservationPoints, see issue 101
                Sim_Theta_U(KT, 1:length(monitorDepthSoilMoisture)) = Theta_UUU(monitorDepthSoilMoisture, 1, KT);
                Sim_Temp(KT, 1:length(monitorDepthTemperature)) = TTT(monitorDepthTemperature, KT);
                Sim_hh(KT, 1:length(monitorDepthSoilMoisture)) = hhh(monitorDepthSoilMoisture, KT);
                Sim_qlh(KT, 1:length(monitorDepthSoilMoisture)) = qlh(monitorDepthSoilMoisture, KT);
                Sim_qlt(KT, 1:length(monitorDepthSoilMoisture)) = qlt(monitorDepthSoilMoisture, KT);
                Sim_qla(KT, 1:length(monitorDepthSoilMoisture)) = qla(monitorDepthSoilMoisture, KT);
                Sim_qvh(KT, 1:length(monitorDepthSoilMoisture)) = qvh(monitorDepthSoilMoisture, KT);
                Sim_qvt(KT, 1:length(monitorDepthSoilMoisture)) = qvt(monitorDepthSoilMoisture, KT);
                Sim_qva(KT, 1:length(monitorDepthSoilMoisture)) = qva(monitorDepthSoilMoisture, KT);
                Sim_qtot(KT, 1:length(monitorDepthSoilMoisture)) = qtot(monitorDepthSoilMoisture, KT);
            end
            if (TEND - TIME) < 1E-3
                for i = 1:ModelSettings.NN
                    hOLD(i) = SoilVariables.h(i);
                    SoilVariables.h(i) = SoilVariables.hh(i);
                    if ModelSettings.Thmrlefc == 1
                        TOLD(i) = SoilVariables.T(i);
                        SoilVariables.T(i) = SoilVariables.TT(i);
                        TTT(i, KT) = SoilVariables.TT(i);
                        TOLD_CRIT(i) = T_CRIT(i);
                        T_CRIT(i) = TT_CRIT(i);
                        hOLD_frez(i) = SoilVariables.h_frez(i);
                        SoilVariables.h_frez(i) = SoilVariables.hh_frez(i);
                    end
                    if ModelSettings.Soilairefc == 1
                        P_gOLD(i) = P_g(i);
                        P_g(i) = P_gg(i);
                    end
                end
            end
        end

        % Recharge calculations, added by Mostafa
        if GroundwaterSettings.GroundwaterCoupling == 1 % Groundwater coupling is enabled
            gwfluxes = groundwater.calculateGroundwaterRecharge(EnergyVariables, SoilVariables, KT, gwfluxes, ModelSettings, GroundwaterSettings);
            if GroundwaterSettings.gw_Dep <= 1 % soil is fully saturated
                gwfluxes.recharge = 0;
            end
        else
            gwfluxes.recharge = 0;
            gwfluxes.indxRchrg = NaN;
        end

        % set SoilVariables for the rest of the loop
        h = SoilVariables.h;
        hh = SoilVariables.hh;
        T = SoilVariables.T;
        TT = SoilVariables.TT;
        hh_frez = SoilVariables.hh_frez;
        h_frez = SoilVariables.h_frez;
        kk = k;
        % Open files for writing
        file_ids = structfun(@(x) fopen(x, 'a'), fnames, 'UniformOutput', false);
        n_col = io.output_data_binary(file_ids, k, xyt, rad, canopy, ScopeParameters, vi, vmax, options, fluxes, ...
                                      meteo, iter, thermal, spectral, gap, profiles, Sim_Theta_U, Sim_Temp, Trap, ...
                                      Evap, WaterStressFactor, WaterPotential, Sim_hh, Sim_qlh, Sim_qlt, Sim_qvh, ...
                                      Sim_qvt, Sim_qla, Sim_qva, Sim_qtot, ForcingData, RS, RWUs, RWUg);
        fclose("all");
    end
end

if strcmp(bmiMode, 'initialize') || strcmp(bmiMode, 'update')
    % Save the required variables to the model state file.
    % NOTE: bmiVarNames are defined in STEMMUS_SCOPE_exe.m
    % save & load is incompatible with Octave. Octave only supports full run mode, not BMI mode.
    save([Output_dir, 'STEMMUS_SCOPE_state.mat'], bmiVarNames{:}, "-v7.3", "-nocompression");
end

if strcmp(bmiMode, 'finalize')
    % Load the workspace to be able to finalize the model.
    load([OutputPath, 'STEMMUS_SCOPE_state.mat']);
end

if strcmp(bmiMode, 'finalize') || strcmp(runMode, 'full')
    disp('Finalizing STEMMUS_SCOPE');
    if options.verify
        io.output_verification(Output_dir);
    end

    io.bin_to_csv(fnames, n_col, k, options, SoilLayer, GroundwaterSettings, FullCSVfiles);
    save([Output_dir, 'output.mat']);
end

% Calculate the total simulatiom time, added by mostafa
end_time = clock;
simtime_min = etime(end_time, start_time) / 60;
simtime_hr = simtime_min / 60;
disp(['Simulation time is : ' num2str(simtime_hr) ' hrs (' num2str(simtime_min) ' minutes)']);
