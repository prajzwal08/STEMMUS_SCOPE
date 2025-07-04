function [resist_out] = resistancesparallel(resist_in)
    %{
    This function calculates the aerodynamic and boundary layer resistances for the soil and vegetation. 
    The vegetation can be modelled as big leaf or as multilayered parallel structure. 
    The overall idea is adopted from the Wallace and Verhoef (2000) model:
    'Modelling interactions in mixed-plant communities: light, water and carbon dioxide', 
    in: Bruce Marshall, Jeremy A. Roberts (ed), 
    'Leaf Development and Canopy Growth', Sheffield Academic Press, UK. ISBN 0849397691
    but adopted to the multilayered parallel structure of the STEMMUS-SCOPE model.
    %}

        % load Constants
    Constants = io.define_constants();
    kappa     = Constants.kappa;
    
    %% parameters
    Cd        = resist_in.Cd;            % Drag coefficient
    u         = resist_in.u;             % Wind speed at reference height
    L         = resist_in.L;             % Monin-Obukhov length
    LAI       = resist_in.LAI;           % Leaf Area Index
    z0m       = resist_in.zo;            % Roughness length for momentum
    z0ms      = 0.01;                    % Soil roughness height [m]
    d         = resist_in.d;             % Zero-plane displacement height
    z         = resist_in.z;             % Reference height
    h         = resist_in.hc;            % Canopy height
    w         = resist_in.w;             % Leaf width
    nl        = resist_in.nl;            % Number of canopy layers
    p         = resist_in.p * 0.1;       % Ambient pressure [hPa --> kPa]
    Ts        = resist_in.Ts;            % Soil surface temperatures [°C]
    Tsh       = Ts(1) + 273.15;          % Shaded soil temperature [K]
    Tsu       = Ts(2) + 273.15;          % Sunlit soil temperature [K]
    po        = 101.3;                   % Standard atmospheric pressure [kPa]
    To        = 273.15;                  % Reference temperature [K]
    kmu0      = 14.6e-6;         %At To and Po, Kinematic viscosity (m2/s).
    % Derived parameters
    zr          = 2.5 * h;                          % Roughness sublayer thickness [m]

    % Layer discretization for multilayer canopy structure
    dh        = (h / nl) * ones(nl, 1);             % Layer thickness [m], equal spacing across nl layers
    hcum      = flip([0; cumsum(dh)]);              % Cumulative height from top of canopy to ground [m], hcum(1) = h (top), hcum(end) = 0 (ground)
    h_layer     = (hcum(1:end-1) + hcum(2:end)) / 2;  % Midpoint height of each layer [m]
    dLAI      = (LAI / nl) * ones(nl, 1);           % Leaf Area Index per layer (uniform distribution)
    dLAIcum   = cumsum(dLAI);                       % Cumulative LAI from top to bottom

    % Wind extinction coefficient (dimensionless), Wallace & Verhoef (2000), Eq. 33
    n       = Cd * LAI / (2 * kappa^2);                % Wind extinction coefficient for entire canopy []
    n_layer  = Cd * dLAIcum ./ (2 * kappa^2);           % Cumulative wind extinction coefficient at each layer []

    %% stability correction for non-neutral conditions
    % neu        = find(L >= -.001 & L <= .001);
    unst        = find(L <  -4);
    st          = find(L >  4E3);

    % stability correction functions, 
    psi_m_z     = psi_m_monin_obukhov(z - d, L, unst, st); %stability correction for momentum at measurement height (z - d)
    psi_h_z     = psi_h_monin_obukhov(z - d, L, unst, st); %stability correction for heat at measurement height (z - d)
    psi_m_h     = psi_m_monin_obukhov(h - d, L, unst, st); %  stability correction for momentum at canopy height (h - d)
    psi_h_zr    = psi_h_monin_obukhov(zr - d, L, unst, st) .* (z >= zr) + psi_h_z .* (z < zr); % Calculate stability correction for heat at reference height (zr - d),
    phi_star_zr = phi_star(zr, zr, d, L, st, unst);% Calculate phi_star correction at reference height (zr)
    phi_star_h  = phi_star(h, zr, d, L, st, unst);% Calculate phi_star correction at canopy height (h)

    ustar       = max(.001, kappa * u ./ (log((z - d) / z0m) - psi_m_z)); % Friction velocity [m s-1] 
    Kh          = kappa * ustar * (zr - d); % Diffusion coefficient for momentum [m2 s-1] 


    %% wind speed at different levels
    u_h   = max(ustar / kappa * (log((h - d) / z0m) - psi_m_h), 0.01);  % Wind speed at canopy height (m/s)
    u_layer = u_h * exp(n * ((h_layer / h) - 1));                          % Wind speed at different layers (m/s)
    u_s = u_h * exp(n * ((z0ms / h) - 1));                          % Wind speed near soil surface at 0.01 m (m/s)
    u_layer   = max(u_layer,0.01);
    u_s  = max(u_s,0.01);

    %% Aerodynamic resistance calculation
    rai = (z > zr) .* (1 ./ (kappa * ustar) .* (log((z - d) / (zr - d))  - psi_h_z   + psi_h_zr)); % Inertial layer above RSL, W&V Eq 41
    rar = 1 ./ (kappa * ustar) .*  ((zr - h) / (zr - d))      - phi_star_zr + phi_star_h; % Roughness sublayer, W&V Eq 39
    rawc = computeCanopyLayerAerodynamicResistance(h, n_layer, Kh, h_layer); % Canopy layer, W&V Eq 42
    rawssh = h * sinh(n) ./ (n * Kh) * (log((exp(n * (h) / h) - 1) / (exp(n * (h) / h) + 1)) - log((exp(n * (z0ms) / h) - 1) / (exp(n * (z0ms) / h) + 1))); %  shaded part of soil.
    rassu = (log(z/z0ms) - psi_m_monin_obukhov(z, L, unst, st) + psi_m_monin_obukhov(z0ms,L,unst,st))^2 / (kappa^2 * u); % sunlit part of soil

    rai  = min(400, rai);
    rar  = min(400, rar);
    rawc = min(400, rawc);
    rawssh = min(400, rawssh);
    rassu = min(400, rassu);

    rac  = rai + rar + rawc; % Total aerodynamic resistance for the plant canopy.
    rassh = rai + rar + rawssh ; % Total aerodynamic resistance for the shaded soil.
    %% BOundary layer resistance
    b1 = 1;
    rbl = b1*70*sqrt(w./u_layer);  %Canopy boundary layer resistance, W&V Eq 31
    
    %For soil
    kmussh  = kmu0 * (po/p) * (Tsh/To)^(1.81); % KInematic viscosity correction for ambient temperature, shaded soil
    kmussu  = kmu0 * (po/p) * (Tsu/To)^(1.81); % KInematic viscosity correction for ambient temperature, sunlit soil

    Ressu   = z0ms*u_s/kmussu;
    Ressh   = z0ms*u_s/kmussh;

    Stssh   = 2.46 * ((Ressh) ^ (1/4)) - log(7.4); % Stanton Number (KB-1)
    Stssu   = 2.46 * ((Ressu) ^ (1/4)) - log(7.4); % Stanton Number (KB-1)

    rbssh = Stssh / (kappa * ustar); % Boundary layer resistance of soil for shaded soil.
    rbssu = Stssu / (kappa * ustar); % Boundary layer resistance of soil for sunlit soil.
  
    % % Output
    % resist_out.rac  = min(4E2, rac);         % to prevent unrealistically high resistances
    % resist_out.rassh = min(4E2, rassh);        % to prevent unrealistically high resistances
    % resist_out.rassu = min(4E2, rassu);        % to prevent unrealistically high resistances

    resist_out.rai = rai;
    resist_out.rar = rar;
    resist_out.rawc = rawc;
    resist_out.rawssh = rawssh;
    resist_out.rac = rac;
    resist_out.rassh = rassh;
    resist_out.rassu = rassu;

    resist_out.rbl = rbl;
    
    resist_out.rbssh = rbssh;% Boundary layer resistance of soil for shaded soil.
    resist_out.rbssu = rbssu; % Boundary layer resistance of soil for sunlit soil.
    
    resist_out.ustar = ustar;
    resist_out.u_h = u_h;
    resist_out.u_s = u_s;
    resist_out.L = L;

    return

%% Stability correction function for momentum (Paulson, 1970)
function pm = psi_m_monin_obukhov(z, L, unst, st)
    % Inputs:
    %   z    - height above displacement height [m]
    %   L    - Monin-Obukhov length [m]
    %   unst - logical indices for unstable conditions (L < 0)
    %   st   - logical indices for stable conditions (L > 0)
    %
    % Output:
    %   pm   - stability correction for momentum

    pm = zeros(size(L));                  % Initialize output array

    % Unstable conditions (L < 0)
    x = (1 - 16 * z ./ L(unst)).^(1/4);
    pm(unst) = 2 * log((1 + x) ./ 2) + log((1 + x.^2) ./ 2) - 2 * atan(x) + pi/2;

    % Stable conditions (L > 0)
    pm(st) = -5 * z ./ L(st);
    return

%% Stability correction function for heat (Paulson, 1970)
function ph = psi_h_monin_obukhov(z, L, unst, st)
    % Inputs:
    %   z    - height above displacement height [m]
    %   L    - Monin-Obukhov length [m]
    %   unst - logical indices for unstable conditions (L < 0)
    %   st   - logical indices for stable conditions (L > 0)
    %
    % Output:
    %   ph   - stability correction for heat

    ph = zeros(size(L));                        % Initialize output array

    % Unstable conditions (L < 0)
    x = (1 - 16 * z ./ L(unst)).^(1/4);
    ph(unst) = 2 * log((1 + x.^2) ./ 2);

    % Stable conditions (L > 0)
    ph(st) = -5 * z ./ L(st);
    return

%% Stability correction function phi_star (Paulson, 1970)
function phs = phi_star(z, zR, d, L, st, unst)
    % Inputs:
    %   z    - measurement height [m]
    %   zR   - reference height [m]
    %   d    - zero-plane displacement height [m]
    %   L    - Monin-Obukhov length [m]
    %   st   - logical indices for stable conditions (L > 0)
    %   unst - logical indices for unstable conditions (L < 0)
    %
    % Output:
    %   phs  - stability correction factor phi_star

    phs = zeros(size(L));                  % Initialize output array

    % Unstable conditions (L < 0)
    x = (1 - 16 * z ./ L(unst)).^0.25;
    phs(unst) = ((z - d) ./ (zR - d)) .* (x.^2 - 1) ./ (x.^2 + 1);

    % Stable conditions (L > 0)
    phs(st) = -5 * z ./ L(st);
    return

function rac = computeCanopyLayerAerodynamicResistance(h, nLayer, Kh, hMid)
    % computeCanopyLayerAerodynamicResistance Calculate aerodynamic resistance per canopy layer using Wallace & Verhoef Eq. 42
    %
    % Inputs:
    %   h       - Canopy height (scalar)
    %   nLayer  - Dimensionless shape parameter (scalar or array)
    %   Kh      - Eddy diffusivity (scalar or array)
    %   hMid    - Midpoint heights of canopy layers (same size as nLayer)
    %
    % Output:
    %   rac     - Aerodynamic resistance for each canopy layer (same size as nLayer)
    %
    % Reference:
    %   Wallace & Verhoef (2000), Eq. 42

    % Compute hyperbolic sine of nLayer
    sinh_nLayer = sinh(nLayer);

    % Calculate numerator and denominator for the expression
    numerator = sinh_nLayer;
    denominator = nLayer .* Kh;

    % Term1: log((exp(nLayer) - 1) / (exp(nLayer) + 1))
    term1 = log((exp(nLayer) - 1) ./ (exp(nLayer) + 1));

    % Calculate zeta = nLayer * (hMid / h)
    zeta = nLayer .* (hMid ./ h);

    % Term2: log((exp(zeta) - 1) / (exp(zeta) + 1))
    term2 = log((exp(zeta) - 1) ./ (exp(zeta) + 1));

    % Calculate aerodynamic resistance rac for each canopy layer
    rac = h .* (numerator ./ denominator) .* (term1 - term2);
    return 