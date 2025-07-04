% this function is edited for multi layer resistance by prajzwal on APril 11, 2025.
function [iter, fluxes, rad, thermal, profiles, soil, RWU, frac, WaterStressFactor,WaterPotential]             ...
         = ebal(iter, options, spectral, rad, gap, leafopt,  ...
                angles, meteo, soil, canopy, leafbio, xyt, k, profiles, Delt_t, ...
                Rl, SoilVariables, VanGenuchten, InitialValues, ModelSettings, GroundwaterSettings,WaterPotential,KT)

    %{
        function ebal.m calculates the energy balance of a vegetated surface

     authors:      Christiaan van der Tol (tol@itc.nl)
                   Joris Timmermans (j_timmermans@itc.nl)
     date          26 Nov 2007 (CvdT)
     updates       29 Jan 2008 (JT & CvdT)     converted into a function
                   11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
                   14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
                   31 Jul 2008 (CvdT)          Included Pntot in output
                   19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
                   07 Nov 2008 (CvdT)          Changed layout
                   18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
                      Feb 2012 (WV)            introduced structures for variables
                      Sep 2013 (JV, CvT)       introduced additional biochemical model

     parent: master.m (script)
     uses:
           RTMt_sb.m, RTMt_planck.m (optional), RTMf.m (optional)
           resistances.m
           heatfluxes.m
           biochemical.m
           soil_respiration.m

     Table of contents of the function

       1. Initialisations for the iteration loop
               intial values are attributed to variables
       2. Energy balance iteration loop
               iteration between thermal RTM and surface fluxes
       3. Write warnings whenever the energy balance did not close
       4. Calculate vertical profiles (optional)
       5. Calculate spectrally integrated energy, water and CO2 fluxes

     The energy balance iteration loop works as follows:

     RTMo              More or less the classic SAIL model for Radiative
                       Transfer of sun and sky light (no emission by the vegetation)
     While continue    Here an iteration loop starts to close the energy
                       balance, i.e. to match the micro-meteorological model
                       and the radiative transfer model
       RTMt_sb         A numerical Radiative Transfer Model for thermal
                       radiation emitted by the vegetation
       resistances     Calculates aerodynamic and boundary layer resistances
                       of vegetation and soil (the micro-meteorological model)
       biochemical     Calculates photosynthesis, fluorescence and stomatal
                       resistance of leaves (or biochemical_MD12: alternative)
       heatfluxes      Calculates sensible and latent heat flux of soil and
                       vegetation
                       Next soil heat flux is calculated, the energy balance
                       is evaluated, and soil and leaf temperatures adjusted
                       to force energy balance closure
     end {while continue}

     meanleaf          Integrates the fluxes over all leaf inclinations
                       azimuth angles and layers, integrates over the spectrum

     usage:
     [iter,fluxes,rad,profiles,thermal]             ...
             = ebal(iter,options,spectral,rad,gap,leafopt,  ...
                    angles,meteo,soil,canopy,leafbio)

     The input and output are structures. These structures are further
     specified in a readme file.

     Input:

       iter        numerical parameters used in the iteration for energy balance closure
       options     calculation options
       spectral    spectral resolutions and wavelengths
       rad         incident radiation
       gap         probabilities of direct light penetration and viewing
       leafopt     leaf optical properties
       angles      viewing and observation angles
       soil        soil properties
       canopy      canopy properties
       leafbio     leaf biochemical parameters

     Output:

       iter        numerical parameters used in the iteration for energy balance closure
       fluxes      energy balance, turbulent, and CO2 fluxes
       rad         radiation spectra
       profiles    vertical profiles of fluxes
       thermal     temperatures, aerodynamic resistances and friction velocity
       sfactor     soil water stress factor
       PSI         leaf water potential
    %}

    %% 1. initialisations and other preparations for the iteration loop
    counter         = 0;              %           Iteration counter of ebal
    maxit           = iter.maxit;
    maxEBer         = iter.maxEBer;
    Wc              = iter.Wc;

    CONT            = 1;              %           is 0 when the calculation has finished
    es_fun      = @(T)6.107 * 10.^(7.5 .* T ./ (237.3 + T)); % saturation vapor pressure
    s_fun       = @(es, T) es * 2.3026 * 7.5 * 237.3 ./ (237.3 + T).^2; % slope of saturation vapor pressure
    t               = xyt.t(k);
    Ta              = meteo.Ta;
    ea              = meteo.ea;
    Ca              = meteo.Ca;
    Ts              = soil.Ts;
    p               = meteo.p;
    if options.soil_heat_method < 2 && options.simulation == 1
        if k > 1                                % 
            Deltat          = Delt_t;           %           Duration of the time interval (s)
        else
            Deltat          = Delt_t;
        end
        x       = [1:12; 1:12]' .* Deltat; %% ?
        Tsold = soil.Tsold; %% Need to check this is this. 
    end

    nl = canopy.nlayers;

    Rnuc  = rad.Rnuc;  % Radiation on sunlit leaves
    GAM   = soil.GAM; % what is this? (Taken from input)
    Tch   = (Ta + .1) * ones(nl, 1);       %           Leaf temperature (shaded leaves) % why random .1 added?
    Tcu   = (Ta + .3) * ones(size(Rnuc)); %           Leaf tempeFrature (sunlit leaves) % why random .3 added?
    ech   = ea * ones(nl, 1);            %           Leaf H2O (shaded leaves)
    ecu   = ea * ones(size(Rnuc));      %           Leaf H2O (sunlit leaves)
    Cch   = Ca * ones(nl, 1);            %           Leaf CO2 (shaded leaves)
    Ccu   = Ca * ones(size(Rnuc));      %           Leaf CO2 (sunlit leaves)
    % Tsold = Ts;                       %           Soil temperature of the previous time step
    L     = -1;                       %           Monin-Obukhov length (Is fixed???)) 
    % load Constants
    Constants = io.define_constants();

    MH2O  = Constants.MH2O;
    Mair  = Constants.Mair;
    rhoa  = Constants.rhoa;
    cp    = Constants.cp;
    g     = Constants.g / 100; % [m s-2] Gravity acceleration
    kappa = Constants.kappa;
    sigmaSB = Constants.sigmaSB;
    Ps    = gap.Ps; % output from RTMo
    nl    = canopy.nlayers; % why again?

    SoilHeatMethod = options.soil_heat_method;
    if ~(options.simulation == 1)
        SoilHeatMethod = 2;
    end

    kV   = canopy.kV;  %% what is kV? %% kV is extinction coefficient for Vcmax in the vertical (maximum at the top). 0 for uniform Vcmax
    xl   = canopy.xl; % division into layer

    % other preparations
    e_to_q          = MH2O / Mair ./ p;  % This is mixing ratio, when it multiplies e, we get vapor pressure.           %           Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
    Fs              = [1 - Ps(end), Ps(end)];      %           Matrix containing values for 1-Ps and Ps of soil
    Fc              = (1 - Ps(1:end - 1))' / nl;      %           Matrix containing values for Ps of canopy, why divide by nl?

    if ~exist('SMCsf', 'var')
        SMCsf = 1;
    end    % HERE COULD BE A STRESS FACTOR FOR VCMAX AS A FUNCTION OF SMC DEFINED
    % but this is at present not
    % incorporated

    fVh             = exp(kV * xl(1:end - 1));
    fVu             = ones(13, 36, nl);
    % Plot

    for i = 1:nl
        fVu(:, :, i) = fVh(i);
    end

    LAI = canopy.LAI;

    %  Get the last time step value (if available)
    % if isfield(WaterPotential, 'leaf') && isempty(WaterPotential.leaf)
    %    PSI = 0;  % Default if no previous step
    % else
    %     PSI = WaterPotential.leaf(end);  % Access last value
    % end

    PSI = 0;

    [bbx] = Max_Rootdepth(InitialValues.bbx, ModelSettings);
                % Reverse BBX
  
    [PSIs, rsss, rrr, rxx] = calc_rsoil(Rl, ModelSettings, SoilVariables, VanGenuchten, bbx, GroundwaterSettings);
    [sfactor] = calc_sfactor(Rl, VanGenuchten.Theta_s, VanGenuchten.Theta_r, SoilVariables.Theta_LL, bbx, Ta, VanGenuchten.Theta_f);
    PSIss = PSIs(ModelSettings.NL, 1);
    %% 2. Energy balance iteration loop

    % 'Energy balance loop (Energy balance and radiative transfer)

    while CONT                          % while energy balance does not close
        
        % 2.1. Net radiation of the components
        % Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
        rad  = RTMt_sb(spectral, rad, soil, leafopt, canopy, gap, angles, Tcu, Tch, Ts(2), Ts(1), 1);
        % Add net radiation of (1) solar and sky and (2) thermal emission model
        
        Rnhct = rad.Rnhct;
        Rnuct = rad.Rnuct;
        Rnhst = rad.Rnhst;
        Rnust = rad.Rnust;

        Rnhc = rad.Rnhc;
        Rnuc = rad.Rnuc;
        Rnhs = rad.Rnhs;
        Rnus = rad.Rnus;

        Rnch        = Rnhc + Rnhct;             %           Canopy (shaded) net radiation
        Rncu        = Rnuc + Rnuct;             %           Canopy (sunlit) net radiation
        Rnsh        = Rnhs + Rnhst;             %           Soil   (shaded) net radiation
        Rnsu        = Rnus + Rnust;             %           Soil   (sunlit) net radiation
        Rns         = [Rnsh Rnsu]';             %           Soil   (sun+sh) net radiation

        % 2.2. Aerodynamic roughness
        % calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
       
        resist_in.u   = max(meteo.u, .2);
        resist_in.L   = L;
        resist_in.LAI = canopy.LAI;
        resist_in.rbs = soil.rbs;
        resist_in.rss = soil.rss;
        resist_in.rwc = canopy.rwc;
        resist_in.zo  = canopy.zo;
        resist_in.d   = canopy.d;
        resist_in.z   = meteo.z;
        resist_in.hc  = canopy.hc;
        resist_in.w   = canopy.leafwidth;
        resist_in.Cd  = canopy.Cd;
        resist_in.nl = canopy.nlayers;
        resist_in.p = meteo.p; % in Kpa, for calculating kinematic viscosity of air.
        resist_in.Ta = meteo.Ta; % In Kelivin , for calculating kinematic viscosity of air.
        resist_in.Ts = Ts; %Soil temperature
        resist_in.Tch = Tch;  % Leaf temperature (shaded)
        resist_in.Tcu = Tcu; % leaf temperature (sunlit)
       
        [resist_out]  = resistancesparallel(resist_in);
       
        
        % 2.3. Biochemical processes

        % photosynthesis (A), fluorescence factor (F), and stomatal resistance (rcw), for shaded (1) and sunlit (h) leaves
        biochem_in.Fluorescence_model = options.Fluorescence_model;
        biochem_in.Type         = leafbio.Type;
        biochem_in.p            = p;
        biochem_in.m            = leafbio.m;
        biochem_in.BallBerry0   = leafbio.BallBerry0;
        biochem_in.O            = meteo.Oa;
        biochem_in.Rdparam      = leafbio.Rdparam;

        if options.Fluorescence_model == 2    % specific for the v.Caemmerer-Magnani model
            biochem_in.Tyear        = leafbio.Tyear;
            biochem_in.beta         = leafbio.beta;
            biochem_in.qLs          = leafbio.qLs;
            biochem_in.NPQs        = leafbio.kNPQs;
            biochem_in.stressfactor = leafbio.stressfactor;
        else
            % specific for Berry-v.d.Tol model
            biochem_in.tempcor      = options.apply_T_corr;
            biochem_in.Tparams      = leafbio.Tparam;
            biochem_in.stressfactor = SMCsf;
        end

        % for shaded leaves
        biochem_in.T        = Tch;
        biochem_in.eb       = ech;
        biochem_in.Vcmo     = fVh .* leafbio.Vcmo;
        biochem_in.Cs       = Cch;
        biochem_in.Q        = rad.Pnh_Cab * 1E6;

        if options.Fluorescence_model == 2
            biochem_out = biochemical_MD12(biochem_in);
        else
            Ci_input = [];
            biochem_out = biochemical(biochem_in, sfactor, Ci_input);
        end

        Ah                  = biochem_out.A;
        Ahh                 = biochem_out.Ag;
        Cih                 = biochem_out.Ci;
        Fh                  = biochem_out.eta;
        rcwh                = biochem_out.rcw;
        qEh                 = biochem_out.qE; % vCaemmerer- Magnani does not generate this parameter (dummy value)
        Knh                 = biochem_out.Kn;

        % for sunlit leaves
        biochem_in.T        = Tcu;
        biochem_in.eb       = ecu;
        biochem_in.Vcmo     = fVu .* leafbio.Vcmo;
        biochem_in.Cs       = Ccu;
        biochem_in.Q        = rad.Pnu_Cab * 1E6;

        if options.Fluorescence_model == 2
            biochem_out = biochemical_MD12(biochem_in);
        else
            Ci_input = [];
            biochem_out = biochemical(biochem_in, sfactor, Ci_input);
        end

        Au                  = biochem_out.A; % Ag? or A?
        Auu                 = biochem_out.Ag;   % GPP calculation.
        Ciu                 = biochem_out.Ci;
        Fu                  = biochem_out.eta;
        rcwu                = biochem_out.rcw;
        qEu                 = biochem_out.qE;
        Knu                 = biochem_out.Kn;

        Pinh                = rad.Pnh;
        Pinu                = rad.Pnu;
        Pinh_Cab            = rad.Pnh_Cab;
        Pinu_Cab            = rad.Pnu_Cab;
        Rnh_PAR             = rad.Rnh_PAR;
        Rnu_PAR             = rad.Rnu_PAR;
        
        % 2.4. Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
        % in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]

        % Resistances
        ustar = resist_out.ustar;
         % Sum of aerodynamic and boundary layer resistance canopy (shaded)
        rach = resist_out.rac + resist_out.rbl; % Sum of aerodynamic and boundary layer resistance canopy (shaded)
        racu  = permute(repmat(resist_out.rac + resist_out.rbl, [1, size(rcwu, 2), size(rcwu, 1)]),[3,2,1]) ; %same shape as rcwu 30*36*13
       
        rss   = soil.rss;  % Stomatal resistance of soil
        % Sum of aerodynamic and boundary layer resistance soil
        ras = [resist_out.rassh + resist_out.rbssh ; resist_out.rassu + resist_out.rbssu];
% Since in Fs, [shaded , sunlight]
        % rac     = (LAI + 1) * (raa + rawc); % why?
        % ras     = (LAI + 1) * (raa + raws); % why? 
        % rach     = raa + rawc; % why?
        % ras     = raws; % why? 
       
        
        for i = 1:30
            
            [lEch, Hch, ech, Cch, lambdah, sh]     = heatfluxes(rach, rcwh, Tch, ea, Ta, e_to_q, PSI, Ca, Cih, es_fun, s_fun);
            [lEcu, Hcu, ecu, Ccu, lambdau, su]     = heatfluxes(racu, rcwu, Tcu, ea, Ta, e_to_q, PSI, Ca, Ciu, es_fun, s_fun);
            [lEs, Hs, ~, ~, lambdas, ss]           = heatfluxes(ras, rss, Ts, ea, Ta, e_to_q, PSIss, Ca, Ca, es_fun, s_fun);
           

            % if any( ~isreal( Cch )) || any( ~isreal( Ccu(:) ))
            %  error('Heatfluxes produced complex values for CO2 concentration!')
            % end

            %  if any( Cch < 0 ) || any( Ccu(:) < 0 )
            %     error('Heatfluxes produced negative values for CO2 concentration!')
            % end
          
            % Hsltot = Fs*Hs;
            % Fc_new = (1 - Ps(1:end - 1)); % fraction of shaded
            % Hcltot = LAI/nl .* (Fc_new.*Hch + (1-Fc_new).*equations.meanleaf(canopy, Hcu, 'angles', Ps)); % This is equivalent to Hctot. 
            
            % Matrix containing values for Ps of canopy, why divide by nl?
            % integration over the layers and sunlit and shaded fractions
            
            Hstot       = Fs * Hs;
            Hctot       = LAI * (Fc * Hch + equations.meanleaf(canopy, Hcu, 'angles_and_layers', Ps));
            
            Htot        = Hstot + Hctot;
            
            %%%%%% Leaf water potential calculate
            lambda1      = (2.501 - 0.002361 * Ta) * 1E6; % This latent heat constant 
            lEctot     = LAI * (Fc * lEch + equations.meanleaf(canopy, lEcu, 'angles_and_layers', Ps)); % latent heat leaves
            if isreal(lEctot) && lEctot < 1000 && lEctot > -300
            else
                lEctot = 0;
            end
            Trans = lEctot / lambda1 / 1000;  % unit: m s-1
            AA1 = PSIs ./ (rsss + rrr + rxx);
            AA2 = 1 ./ (rsss + rrr + rxx);
            BB1 = AA1(~isnan(AA1));
            BB2 = AA2(~isinf(AA2));
            PSI1 = (sum(BB1) - Trans) / sum(BB2); % updated 

            % if isnan(PSI1)
            %     warning('PSI is NaN!')
            %      % Create a dynamic filename based on the value of k
            %     filename = sprintf('/Users/prajzwal/PhD/debugtrial/debugtrial4_%d.mat', k);
            %     SMC_trial = SoilVariables.Theta_LL(1:end - 1, 2);
            %     save(filename, 'PSI1', 'PSIs', 'rsss', 'rrr', 'rxx', 'SMC_trial'); 
            %     PSI1 = -1;
            % end
        
            % if ~isreal(PSI1)
            %     warning('PSI is a complex number!')
            %     filename = sprintf('/Users/prajzwal/PhD/debugtrial/debugtrial4_%d.mat', k);
            %     SMC_trial = SoilVariables.Theta_LL(1:end - 1, 2);
            %     save(filename, 'PSI1', 'PSIs', 'rsss', 'rrr', 'rxx', 'SMC_trial'); 
            %     PSI1 = -1;
            % end

            if abs(PSI - PSI1) < 0.01
                break

            end
            PSI  = (PSI + PSI1) / 2;
        end
        
       
        %%%%%%%
        if SoilHeatMethod == 2
            G = 0.30 * Rns;
        else
            G = GAM / sqrt(pi) * 2 * sum(([Ts'; Tsold(1:end - 1, :)] - Tsold) ./ Deltat .* (sqrt(x) - sqrt(x - Deltat)));
            G = G';
        end
        % 2.5. Monin-Obukhov length L
        L           = -rhoa * cp * ustar.^3 .* (Ta + 273.15) ./ (kappa * g * Htot);           % [1]
        L(L < -1E3)   = -1E3;                                                     % [1]
        L(L > 1E2)    =  1E2;                                                     % [1]
        L(isnan(L)) = -1;                                                       % [1]
       
        % 2.6. energy balance errors, continue criterion and iteration counter
        EBerch      = Rnch - lEch - Hch;
        EBercu      = Rncu - lEcu - Hcu;
        EBers       = Rns  - lEs  - Hs - G;

        counter     = counter + 1;                   %        Number of iterations
        maxEBercu   = max(max(max(abs(EBercu))));
        maxEBerch   = max(abs(EBerch));
        maxEBers    = max(abs(EBers));

        CONT        = (maxEBercu >   maxEBer    | ...
                       maxEBerch >   maxEBer    | ...
                       maxEBers  >   maxEBer)    & ...
                        counter   <   maxit + 1; %        Continue iteration?
        if counter == 5
            Wc = 0.6;
        end
        if counter == 10
            Wc = 0.4;
        end
        if counter == 30
            Wc = 0.1;
        end
       
       
        % 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
        Tch         = Tch + Wc * EBerch ./ ((rhoa * cp) ./ rach + rhoa * lambdah * e_to_q .* sh ./ (rach + rcwh) + 4 * leafbio.emis * sigmaSB * (Tch + 273.15).^3);
        Tcu         = Tcu + Wc * EBercu ./ ((rhoa * cp) ./ racu + rhoa * lambdau * e_to_q .* su ./ (racu + rcwu) + 4 * leafbio.emis * sigmaSB * (Tcu + 273.15).^3);
        Ts          = Ts + Wc * EBers ./ (rhoa * cp ./ ras + rhoa * lambdas * e_to_q .* ss ./ (ras + rss) + 4 * (1 - soil.rs_thermal) * sigmaSB * (Ts + 273.15).^3);
        if ~isreal(Ts)
                keyboard;
        end 
        % Ts contains shaded soil temperature and sunlit soil temperature
        Tch(abs(Tch) > 100) = Ta;
        Tcu(abs(Tcu) > 100) = Ta;
        Ts(abs(Ts) > 100) = Ta;
        if any(isnan(Tch)) || any(isnan(Tcu(:)))
            warning('Canopy temperature gives NaNs');
        end
        if any(isnan(Ts))
            warning('Soil temperature gives NaNs');
        end

    end
   

    iter.counter = counter;
    profiles.etah = Fh;
    profiles.etau = Fu;
    
   
    if SoilHeatMethod < 2
        Tsold(2:end, :) = soil.Tsold(1:end - 1, :);
        Tsold(1, :)  = Ts(:);
        if isnan(Ts)
            Tsold(1, :) = Tsold(2, :);
        end
        if isreal(Ts)
            soil.Tsold = Tsold;
        else
            Tsold(1, :) = Tsold(2, :);
            soil.Tsold = Tsold;
        end
    end

    Tbr         = (rad.Eoutte / Constants.sigmaSB)^0.25;
    Lot_        = equations.Planck(spectral.wlS', Tbr);
    rad.LotBB_  = Lot_;           % Note that this is the blackbody radiance!

    %% 3. Print warnings whenever the energy balance could not be solved
    if counter >= maxit
        warning(['\n warning: maximum number of iteratations exceeded', ...
                 '\n Energy balance error sunlit vegetation = %4.2f W m-2 ', ...
                 '\n Energy balance error shaded vegetation = %4.2f W m-2 ', ...
                 '\n Energy balance error soil              = %4.2f W m-2 '], maxEBercu, maxEBerch, maxEBers);


    end

    %% 4. Calculate the output per layer
    if options.calc_vert_profiles
        [Hcu1d]           = equations.meanleaf(canopy, Hcu,          'angles');   % [nli,nlo,nl]      mean sens heat sunlit leaves
        [lEcu1d]           = equations.meanleaf(canopy, lEcu,         'angles');   % [nli,nlo,nl]      mean latent sunlit leaves
        [Au1d]           = equations.meanleaf(canopy, Au,           'angles');   % [nli,nlo,nl]      mean phots sunlit leaves
        [Fu_Pn1d]           = equations.meanleaf(canopy, Fu .* Pinu_Cab, 'angles');   % [nli,nlo,nl]      mean fluor sunlit leaves
        [qEuL]           = equations.meanleaf(canopy, qEu,          'angles');   % [nli,nlo,nl]      mean fluor sunlit leaves
        % [Pnu1d  ]           = equations.meanleaf(canopy,Pinu,         'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
        % [Pnu1d_Cab  ]       = equations.meanleaf(canopy,Pinu_Cab,     'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
        [Rnu1d]           = equations.meanleaf(canopy, Rncu,         'angles');   % [nli,nlo,nl]      mean net PAR sunlit leaves
        [Tcu1d]           = equations.meanleaf(canopy, Tcu,          'angles');   % [nli,nlo,nl]      mean temp sunlit leaves

        profiles.Tchave     = mean(Tch);                                           % [1]               mean temp shaded leaves
        profiles.Tch        = Tch;                                                 % [nl]
        profiles.Tcu1d      = Tcu1d;                                               % [nl]
        profiles.Tc1d       = (1 - Ps(1:nl)) .* Tch       + Ps(1:nl) .* (Tcu1d);         % [nl]              mean temp leaves, per layer
        profiles.Hc1d       = (1 - Ps(1:nl)) .* Hch       + Ps(1:nl) .* (Hcu1d);         % [nl]              mean sens heat leaves, per layer
        profiles.lEc1d      = (1 - Ps(1:nl)) .* lEch      + Ps(1:nl) .* (lEcu1d);        % [nl]              mean latent heat leaves, per layer
        profiles.A1d        = (1 - Ps(1:nl)) .* Ah        + Ps(1:nl) .* (Au1d);          % [nl]              mean photos leaves, per layer
        profiles.F_Pn1d     = ((1 - Ps(1:nl)) .* Fh .* Pinh_Cab + Ps(1:nl) .* (Fu_Pn1d));  % [nl]           mean fluor leaves, per layer
        profiles.qE         = ((1 - Ps(1:nl)) .* qEh      + Ps(1:nl) .* (qEuL));         % [nl]           mean fluor leaves, per layer
        % profiles.Pn1d       = ((1-Ps(1:nl)).*Pinh     + Ps(1:nl).*(Pnu1d));        %[nl]           mean photos leaves, per layer
        % profiles.Pn1d_Cab   = ((1-Ps(1:nl)).*Pinh_Cab + Ps(1:nl).*(Pnu1d_Cab));        %[nl]           mean photos leaves, per layer
        profiles.Rn1d       = ((1 - Ps(1:nl)) .* Rnch     + Ps(1:nl) .* (Rnu1d));        % [nl]
    end

    %% 5. Calculate spectrally integrated energy, water and CO2 fluxes
    % sum of all leaves, and average leaf temperature
    %     (note that averaging temperature is physically not correct...)
    
    Rnctot          = LAI * (Fc * Rnch + equations.meanleaf(canopy, Rncu, 'angles_and_layers', Ps)); % net radiation leaves
    lEctot          = LAI * (Fc * lEch + equations.meanleaf(canopy, lEcu, 'angles_and_layers', Ps)); % latent heat leaves
    Hctot           = LAI * (Fc * Hch  + equations.meanleaf(canopy, Hcu, 'angles_and_layers', Ps)); % sensible heat leaves
    % Actot           = LAI*(Fc*Ah   + equations.meanleaf(canopy,Au  ,'angles_and_layers',Ps)); % photosynthesis leaves
    Actot           = LAI * (Fc * Ahh   + equations.meanleaf(canopy, Auu, 'angles_and_layers', Ps)); % photosynthesis leaves
    Tcave           =  (Fc * Tch  + equations.meanleaf(canopy, Tcu, 'angles_and_layers', Ps)); % mean leaf temperature
    Pntot           = LAI * (Fc * Pinh + equations.meanleaf(canopy, Pinu, 'angles_and_layers', Ps)); % net PAR leaves
    Pntot_Cab       = LAI * (Fc * Pinh_Cab + equations.meanleaf(canopy, Pinu_Cab, 'angles_and_layers', Ps)); % net PAR leaves
    Rntot_PAR       = LAI * (Fc * Rnh_PAR  + equations.meanleaf(canopy, Rnu_PAR, 'angles_and_layers', Ps)); % net PAR leaves
    aPAR_Cab_eta        = LAI * (Fc * (profiles.etah .* Rnh_PAR) + equations.meanleaf(canopy, profiles.etau .* Rnu_PAR, 'angles_and_layers', Ps));
    % ... green ePAR * relative fluorescence emission efficiency
    % [Delta_Rltot] = Root_properties(Rl, Ac, rroot, frac, bbx, KT, DeltZ, sfactor, LAI_msr);
    % Delta_Rl = fc*Delta_Rltot;
    % Rl = Rl + Delta_Rl;
    % Rltot = sum(sum(Rl));
    % fc = Rl./Rltot;
    % sum of soil fluxes and average temperature
    %   (note that averaging temperature is physically not correct...)
    Rnstot          = Fs * Rns;           %                   Net radiation soil
    lEstot          = Fs * lEs;           %                   Latent heat soil
    % Hstot          = Fs*Hs;            %                   Sensible heat soil
    Gtot            = Fs * G;             %                   Soil heat flux
    Tsave           = Fs * Ts;            %                   Soil temperature
    Resp            = Fs * equations.soil_respiration(Ts); %             Soil respiration

    % total fluxes (except sensible heat), all leaves and soil
    Atot            = Actot;            %                   GPP
    Rntot           = Rnctot + Rnstot;  %                   Net radiation
    lEtot           = lEctot + lEstot;  %                   Latent heat
    % Htot           = Hctot  + Hstot;   %                   Sensible heat

    fluxes.Rntot    = Rntot;  % [W m-2]             total net radiation (canopy + soil)
    fluxes.lEtot    = lEtot;  % [W m-2]             total latent heat flux (canopy + soil)
    fluxes.Htot     = Htot;   % [W m-2]             total sensible heat flux (canopy + soil)
    fluxes.Atot     = Atot;   % [umol m-2 s-1]      total net CO2 uptake (canopy + soil)
    fluxes.Rnctot   = Rnctot; % [W m-2]             canopy net radiation
    fluxes.lEctot   = lEctot; % [W m-2]             canopy latent heat flux
    fluxes.Hctot    = Hctot;  % [W m-2]             canopy sensible heat flux
    fluxes.Actot    = Actot;  % [umol m-2 s-1]      canopy net CO2 uptake
    fluxes.Rnstot   = Rnstot; % [W m-2]             soil net radiation
    fluxes.lEstot   = lEstot; % [W m-2]             soil latent heat flux
    fluxes.Hstot    = Hstot;  % [W m-2]             soil sensible heat flux
    fluxes.Gtot     = Gtot;   % [W m-2]             soil heat flux
    if isnan(fluxes.Rntot)
        fluxes.lEtot = Rntot;
        fluxes.Htot = Rntot;
        fluxes.Atot = Rntot;
        fluxes.Rnctot   = Rntot; % [W m-2]             canopy net radiation
        fluxes.lEctot   = Rntot; % [W m-2]             canopy latent heat flux
        fluxes.Hctot    = Rntot;  % [W m-2]             canopy sensible heat flux
        fluxes.Actot    = Rntot;  % [umol m-2 s-1]      canopy net CO2 uptake
        fluxes.Rnstot   = Rntot; % [W m-2]             soil net radiation
        fluxes.lEstot   = Rntot; % [W m-2]             soil latent heat flux
        fluxes.Hstot    = Rntot;  % [W m-2]             soil sensible heat flux
        fluxes.Gtot     = Rntot;   % [W m-2]             soil heat flux
    end
    fluxes.Resp     = Resp;   % [umol m-2 s-1]      soil respiration
    fluxes.aPAR     = Pntot;  % [umol m-2 s-1]      absorbed PAR
    fluxes.aPAR_Cab = Pntot_Cab; % [umol m-2 s-1]      absorbed PAR
    fluxes.aPAR_Wm2 = Rntot_PAR; % [W m-2]      absorbed PAR
    fluxes.aPAR_Cab_eta = aPAR_Cab_eta;
    fluxes.GPP    = Actot * 12 / 1000000000;  % [kg C m-2 s-1]      gross primary production
    fluxes.NEE      = (Resp - Actot) * 12 / 1000000000; % [kg C m-2 s-1]      net primary production

    thermal.Ta    = Ta;       % [oC]                air temperature (as in input)
    thermal.Ts    = Ts;       % [oC]                soil temperature, sunlit and shaded [2x1]
    thermal.Tcave = Tcave;    % [oC]                weighted average canopy temperature
    thermal.Tsave = Tsave;    % [oC]                weighted average soil temperature
    thermal.raa   = resist_out.rar;      % [s m-1]             total aerodynamic resistance above canopy
    thermal.rawc  = mean(rach);     % [s m-1]             aerodynamic resistance below canopy for canopy
    thermal.raws  = resist_out.rawssh;     % [s m-1]             aerodynamic resistance below canopy for soil
    thermal.ustar = ustar;    % [m s-1]             friction velocity
    thermal.Tcu   = Tcu;
    thermal.Tch   = Tch;
    fluxes.Au     = Au;
    fluxes.Ah     = Ah;

    RWU = (PSIs - PSI) ./ (rsss + rrr + rxx) .* bbx;
    nn = numel(RWU);
    for i = 1:nn
        if isnan(RWU(i))
            RWU(i) = 0;
        end
    end
    for i = 1:nn
        if RWU(i) < 0
            RWU(i) = 1 * 1e-20;
        end
    end
    frac = RWU ./ abs(sum(sum(RWU)));
    RWU = (PSIs - PSI) ./ (rsss + rrr + rxx) .* bbx;
    RWU = real(RWU);
    for i = 1:nn
        if isnan(RWU(i))
            RWU(i) = 0;
        end
    end

    profiles.Knu     = Knu;
    profiles.Knh     = Knh;
    % function Tnew = update(Told, Wc, innovation)
    %     Tnew        = Wc.*innovation + (1-Wc).*Told;
    % return
    
    WaterStressFactor.soil = sfactor;
    WaterPotential.leaf = PSI;
    % Open the file for the first time (write mode)
   % Profiles
   
    

    fid = fopen('/home/khanalp/paper01/output/vMLWV/IT-Lav/anciliary_output.csv', 'a');  % 'a' mode to append; if file doesn't exist, it'll create one
    
    if fid == -1
        error('Error opening the file for writing.');
    end
    
    % Write the header only if KT == 1 and the file is empty
    if ftell(fid) == 0  % Check if file is empty by checking the position of the file pointer
        fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
            'KT','u', 'L', 'LAI', 'z', 'h','w', ...
            'ustar','u_h','u_s','d', ...
            'rai','rar','rawssh', ...
            'rassh', 'rassu', ...
            'rbssh','rbssu',...
            'Fsh','Tsh','Tsu','Tsave','Tcave', 'Rnssh','Rnssu','Hssh','Hssu','lEssh','lEssu','Gssh','Gssu');

    end 
    
    fprintf(fid, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', ...
        KT,resist_in.u, resist_out.L, resist_in.LAI, resist_in.z, resist_in.hc, resist_in.w, ...
        resist_out.ustar, resist_out.u_h ,resist_out.u_s,  resist_in.d, ...
        resist_out.rai, resist_out.rar,resist_out.rawssh, ...
        resist_out.rassh, resist_out.rassu, ...
        resist_out.rbssh, resist_out.rbssu, ...
        Fs(1), Ts(1), Ts(2), Tsave, Tcave, Rns(1), Rns(2), Hs(1), Hs(2), lEs(1), lEs(2), G(1), G(2));

    
    % Close the file after the operation
    fclose(fid);
   % Define output file paths
    rac_file         = '/home/khanalp/paper01/output/vMLWV/IT-Lav/rac.csv';
    rbl_file      = '/home/khanalp/paper01/output/vMLWV/IT-Lav/rbl.csv';
    rcwh_file     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/rcwh.csv';
    % rbcsh_file     = '/home/khanalp/paper01/output/vMLWV/IT-MBo/rbcsh.csv';
    fraction_file    = '/home/khanalp/paper01/output/vMLWV/IT-Lav/fraction.csv';
    profile_Hcu      = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_Hcu.csv';
    profile_Rncu     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_Rncu.csv';
    profile_lEcu     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_lEcu.csv';
    profile_Tcu     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_Tcu.csv';
    profile_Hch      = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_Hch.csv';
    profile_Rnch     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_Rnch.csv';
    profile_lEch     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_lEch.csv';
    profile_Tch     = '/home/khanalp/paper01/output/vMLWV/IT-Lav/profile_Tch.csv';


    % Compute derived quantities
    fraction = 1 - Ps(1:end - 1);  % shaded fraction [30x1]
    Hlcu     = equations.meanleaf(canopy, Hcu,  'angles', Ps);   % [30x1]
    Rnlcu    = equations.meanleaf(canopy, Rncu, 'angles', Ps);   % [30x1]Rnlcu    = equations.meanleaf(canopy, Rncu, 'angles', Ps);   % [30x1]
    lElcu    = equations.meanleaf(canopy, lEcu, 'angles', Ps);   % [30x1]
    Tlcu     = equations.meanleaf(canopy, Tcu, 'angles', Ps);   % [30x1]
   
    % Save all variables using the helper function
    append_or_create_csv(rac_file,      resist_out.rac(:));
    append_or_create_csv(rbl_file,   resist_out.rbl(:));
    append_or_create_csv(fraction_file, fraction(:));
    append_or_create_csv(profile_Hcu,   Hlcu(:));
    append_or_create_csv(profile_Rncu,  Rnlcu(:));
    append_or_create_csv(profile_lEcu,  lElcu(:));
    append_or_create_csv(profile_Tcu,   Tlcu(:));
    append_or_create_csv(profile_Hch,   Hch(:));
    append_or_create_csv(profile_Rnch,  Rnch(:));
    append_or_create_csv(profile_lEch,  lEch(:));
    append_or_create_csv(profile_Tch,   Tch(:));
    append_or_create_csv(rcwh_file, rcwh(:));
    % append_or_create_csv(rbcsh_file, resist_out.rbcsh(:));

    % if KT == 25
    %     keyboard;
    % end

% Define utility function to append or create CSV
function append_or_create_csv(file_path, column_data)
    if isfile(file_path)
        existing_data = readmatrix(file_path);
        updated_data = [existing_data, column_data];
    else
        updated_data = column_data;
    end
    writematrix(updated_data, file_path);
return
