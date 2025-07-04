function [resist_out] = resistancesHF(resist_in)
%==========================================================================
% resistancesHF - Compute aerodynamic and boundary layer resistances
%                 for high-frequency land-atmosphere modeling using Harman
%                 and Finigan , 2007 and 2008
%
% INPUT:
%   resist_in - struct containing input parameters:
%       Cd   - drag coefficient
%       u    - wind speed (m/s)
%       L    - Obukhov length (m)
%       LAI  - leaf area index
%       z    - measurement height (m)
%       h    - canopy height (m)
%       w    - leaf width (m)
%       nl   - number of canopy layers
%       p    - pressure (kPa)
%       Ta   - air temperature (K)
%
% OUTPUT:
%   resist_out - struct containing resistance components and ustar
%==========================================================================

% Load constants
Constants = io.define_constants();
kappa = Constants.kappa;

%--------------------------------------------------------------------------
% Unpack input
    
    u         = resist_in.u;             % Wind speed at reference height
    L         = resist_in.L;             % Monin-Obukhov length
    LAI       = resist_in.LAI;           % Leaf Area Index
    z0ms      = 0.01;                    % Soil roughness height [m]
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

    % Constants and empirical coefficients related to Harman and Finigan, 2007
    c2 = 0.5;
    beta_neutral_max = 0.35;

    % Layer discretization for multilayer canopy structure
    dh        = (h / nl) * ones(nl, 1);             % Layer thickness [m], equal spacing across nl layers
    hcum      = flip([0; cumsum(dh)]);              % Cumulative height from top of canopy to ground [m], hcum(1) = h (top), hcum(end) = 0 (ground)
    h_layer     = (hcum(1:end-1) + hcum(2:end)) / 2;  % Midpoint height of each layer [m]
    % dLAI      = (LAI / nl) * ones(nl, 1);           % Leaf Area Index per layer (uniform distribution)
    % % dLAIcum   = cumsum(dLAI);                       % Cumulative LAI from top to bottom

%--------------------------------------------------------------------------
% Load lookup tables for ψ_hat (roughness sublayer corrections)
    lookup_psihat = "psihat.nc";
    dtLgridM = ncread(lookup_psihat, 'dtLgridM');
    zdtgridM = ncread(lookup_psihat, 'zdtgridM');
    psigridM = ncread(lookup_psihat, 'psigridM');
    dtLgridH = ncread(lookup_psihat, 'dtLgridH');
    zdtgridH = ncread(lookup_psihat, 'zdtgridH');
    psigridH = ncread(lookup_psihat, 'psigridH');

%--------------------------------------------------------------------------
    % Calculate canopy length scale
    Lc = 4 * h / LAI;

    % Compute neutral beta and limit it
    c_beta = kappa^2 * (log((h + z0ms) / z0ms))^2;
    beta_neutral = min(sqrt(c_beta + 0.3 * LAI), beta_neutral_max);

    % Avoid near-zero L
    if abs(L) <= 0.1
        L = 0.1;
    end

    LcL = Lc / L;

    % Compute beta = u*/u(h)
    beta = GetBeta(beta_neutral, LcL);

    % Compute Prandtl and Schmidt numbers
    [Pr, Sc] = GetPrSc(beta_neutral, beta_neutral_max, LcL);

%--------------------------------------------------------------------------
    % Displacement height (d) from top of canopy
    dt = beta^2 * Lc * (1 - exp(-0.25 * LAI / beta^2));
    dt = min(h, dt);
    d = h - dt;

    % Limit Obukhov length to avoid extreme cases, the phim, phic, psim, psic
    % are only valid from -2 to 1, Foken, 2006. 
    zeta = (z - d) / L;
    if zeta >= 0
        zeta = min(1.0, max(zeta, 0.01));
    else
        zeta = max(-2.0, min(zeta, -0.01));
    end
    L = (z - d) / zeta;

%--------------------------------------------------------------------------
    % Stability functions
    phim = phi_m_monin_obukhov((h - d) / L);
    % phic = phi_c_monin_obukhov((h - d) / L);

    % Roughness sublayer correction coefficients
    c1_m = (1 - kappa / (2 * beta * phim)) * exp(0.5 * c2); %For momentum
    % c1_c = (1 - (Sc * kappa) / (2 * beta * phic)) * exp(0.5 * c2); %For heat,vapor

% Deviations from log-law due to roughness sublayer
    dv_m = - psi_m_monin_obukhov((z - d)/(L)) + ...
            psi_m_monin_obukhov((h - d)/L) + ...
            c1_m * (LookupPsihat((z - h)/(h - d), (h - d)/L, zdtgridM, dtLgridM, psigridM) - ...
            LookupPsihat((h - h) / (h - d), (h - d) / L, zdtgridM, dtLgridM, psigridM)) + kappa/beta;

    % dv_c = - psi_c_monin_obukhov((z - d)/(h - d)) + ...
    %          psi_c_monin_obukhov((h - d)/L) + ...
    %          c1_c * (LookupPsihat((z - h)/(h - d), (h - d)/L, zdtgridH, dtLgridH, psigridH) - ...
    %          LookupPsihat((h - h) / (h- d), (h- d) / L, zdtgridH, dtLgridH, psigridH));

    % z0m = (h-d) * exp(-k/beta) * exp(-psi_m_monin_obukhov((h-d)/L)+psi_m_monin_obukhov(z))

%--------------------------------------------------------------------------
    % Wind profile and turbulence characteristics
    zlog = log((z - d) / (hc - d));
    ustar = u * kappa / (zlog + dv_m);
    u_h = ustar / beta;
    lm = 2 * beta^3 * Lc; %Mixing length, formula (6, H&F 2007)
    u_layer = u_h * exp((h_layer - h) / (lm / beta)); % Wind speed at each layer
    u_s = u_h * exp((z0ms - h) / (lm / beta)); % %Wind speed at soil height. 
    u_layer   = max(u_layer,0.01);
    u_s  = max(u_s,0.01);

%--------------------------------------------------------------------------
    % Aerodynamic and boundary layer resistances
    rar = (1 / (kappa * ustar)) * (zlog + dv_m); % Resitance from the top of the canopy to the flux tower height.
    rawc = (Sc / (beta * ustar)) * ...
        (exp(-(hmid - h) / (lm / beta)) - exp(-(h- h) / (lm / beta))); % Within canopy resistance, second term is zero, because to canopy height.
    rawssh = (Sc / (beta * ustar)) * ...
        (exp(-(z0ms - h) / (lm / beta)) - exp(-(h- h) / (lm / beta)));% Within canopy resistance for the soil part.
    rassu = (log(z/ z0ms) - psi_m_monin_obukhov(z/L) + ...
            psi_m_monin_obukhov(z0ms/L))/(kappa^2*u); % Aerodynamic resistance from bare soil to reference height

    rac = rar+rawc;
    rassh = rar+rawssh; 
    
    rar  = min(400, rar);
    rawc = min(400, rawc);
    rawssh = min(400, raws);
    rassu = min(400, rassu);

%--------------------------------------------------------------------------
    % Boundary layer part
    rbl = 70 * sqrt(w ./ u_layer);   % Leaf boundary layer resistance
    % rbc = rbl ./ dLAI;          % Leaf to canopy
     %For soil
    kmussh  = kmu0 * (po/p) * (Tsh/To)^(1.81); % KInematic viscosity correction for ambient temperature, shaded soil
    kmussu  = kmu0 * (po/p) * (Tsu/To)^(1.81); % KInematic viscosity correction for ambient temperature, sunlit soil

    %Reynold's number
    Ressu   = z0ms*u_s/kmussu;
    Ressh   = z0ms*u_s/kmussh;

    Stssh   = 2.46 * ((Ressh) ^ (1/4)) - log(7.4); % Stanton Number (KB-1)
    Stssu   = 2.46 * ((Ressu) ^ (1/4)) - log(7.4); % Stanton Number (KB-1)

    rbssh = Stssh / (kappa * ustar); % Boundary layer resistance of soil for shaded soil.
    rbssu = Stssu / (kappa * ustar); % Boundary layer resistance of soil for sunlit soil.


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Pack output
    resist_out.rar = rar;
    resist_out.rawc = rawc;
    resist_out.rawssh = rawssh;
    resist_out.rac = rac;
    resist_out.rassh = rassh;
    resist_out.rassu = rassu;

    resist_out.rbl = rbl;
    resist_out.rbssh = rbssh;
    resist_out.rbssu = rbssu;

    
    resist_out.ustar = ustar;
    resist_out.u_h = u_h;
    resist_out.u_s = u_s;
    resist_out.L = L;

    resist_out.lm = lm;
    resist_out.beta = beta;
    resist_out.d = d;
    resist_out.Lc = Lc;
    return

function beta = GetBeta(beta_neutral, LcL)
%GETBETA Computes beta = u*/u(h) for current scale length (LcL)
%
% Inputs:
%   beta_neutral - Neutral value for beta = u*/u(h), ~0.3-0.5
%   LcL          - Canopy scale Lc divided by Obukhov length
%
% Output:
%   beta         - Value of u*/u(h) at canopy top

    % Copy LcL and use local variable
    LcL_val = LcL;

    if LcL_val <= 0
        % Unstable case: quadratic equation for beta^2
        bb = 16 * LcL_val * beta_neutral^4;
        beta = sqrt(0.5 * (-bb + sqrt(bb^2 + 4 * beta_neutral^4)));
    else
        % Stable case: cubic equation for beta
        aa = 5 * LcL_val;
        bb = 0;
        cc = 1;
        dd = -beta_neutral;

        qq = (2 * bb^3 - 9 * aa * bb * cc + 27 * (aa^2) * dd)^2 ...
           - 4 * (bb^2 - 3 * aa * cc)^3;
        qq = sqrt(qq);
        rr = 0.5 * (qq + 2 * bb^3 - 9 * aa * bb * cc + 27 * (aa^2) * dd);
        rr = rr^(1/3);
        beta = -(bb + rr) / (3 * aa) - (bb^2 - 3 * aa * cc) / (3 * aa * rr);
    end
    
   
    % Limit beta between 0.20 and 0.50
    beta = min(0.50, max(beta, 0.20));
    return

function [Pr,Sc] = GetPrSc(beta_neutral, beta_neutral_max, LcL)
%GETPRSC Computes the Prandtl number (Pr) at canopy top for current Obukhov length (LcL)
%
% Inputs:
%   beta_neutral      - Neutral value for beta = u*/u(h), ~0.3-0.5
%   beta_neutral_max  - Maximum neutral value for beta
%   LcL               - Canopy scale Lc divided by Obukhov length
%
% Output:
%   Pr                - Prandtl number at canopy top

    % Constants
    Prn  = 0.5;   % Neutral value for Pr
    Prvr = 0.3;   % Magnitude of variation of Pr with stability
    Prsc = 2.0;   % Scale of variation of Pr with stability

    Scn  = Prn;   % Neutral value for Sc
    Scvr = Prvr;  % Magnitude of variation of Sc with stability
    Scsc = Prsc;  % Scale of variation of Sc with stability

    % Compute Prandtl number
    Pr = Prn + Prvr * tanh(Prsc * LcL);
    Pr = (1 - beta_neutral / beta_neutral_max) * 1 + (beta_neutral / beta_neutral_max) * Pr;

    % Compute Schmidt number (if needed, assign or return as output)
    Sc = Scn + Scvr * tanh(Scsc * LcL);
    Sc = (1 - beta_neutral / beta_neutral_max) * 1 + (beta_neutral / beta_neutral_max) * Sc;
    return

function psihat = LookupPsihat(zdt, dtL, zdtgrid, dtLgrid, psigrid)
    % LookupPsihat determines the RSL function psihat from a lookup table
    % using bilinear interpolation for input zdt and dtL.
    %
    % Inputs:
    %   zdt      - Height (above canopy) normalized by dt (scalar)
    %   dtL      - dt/L (displacement height / Obukhov length) (scalar)
    %   zdtgrid  - nZ x 1 array of zdt values
    %   dtLgrid  - 1 x nL array of dtL values
    %   psigrid  - nL x nZ array of psihat values
    %STEMMUS_SCOPE_model/STEMMUS_SCOPE_v2/STEMMUS_SCOPE/src/psihat.nc
    % Output:
    %   psihat   - Interpolated psihat value

    nL = length(dtLgrid);
    nZ = length(zdtgrid);


    % Handle edge cases for dtL
    if dtL <= dtLgrid(1)
        L1 = 1; L2 = 1;
        wL1 = 0.5; wL2 = 0.5;
    elseif dtL > dtLgrid(end)
        L1 = nL; L2 = nL;
        wL1 = 0.5; wL2 = 0.5;
    else
        L1 = 0; L2 = 0;
        for jj = 1:(nL-1)
            if dtL > dtLgrid(jj) && dtL <= dtLgrid(jj+1)
                L1 = jj;
                L2 = jj+1;
                denom = dtLgrid(L2) - dtLgrid(L1);
                wL1 = (dtLgrid(L2) - dtL) / denom;
                wL2 = 1.0 - wL1;
                break;
            end
        end
        if L1 == 0 || L2 == 0
            error('ERROR: LookupPsihat - Indices L1 and L2 not found.');
        end
    end

    % Handle edge cases for zdt
    if zdt > zdtgrid(1)
        Z1 = 1; Z2 = 1;
        wZ1 = 0.5; wZ2 = 0.5;
    elseif zdt < zdtgrid(end)
        Z1 = nZ; Z2 = nZ;
        wZ1 = 0.5; wZ2 = 0.5;
    else
        Z1 = 0; Z2 = 0;
        for ii = 1:(nZ-1)
            if zdt < zdtgrid(ii) && zdt >= zdtgrid(ii+1)
                Z1 = ii;
                Z2 = ii+1;
                denom = zdtgrid(Z1) - zdtgrid(Z2);
                wZ1 = (zdt - zdtgrid(Z2)) / denom;
                wZ2 = 1.0 - wZ1;
                break;
            end
        end
        if Z1 == 0 || Z2 == 0
            error('ERROR: LookupPsihat - Indices Z1 and Z2 not found.');
        end
    end

    % fprintf('L1 = %d, L2 = %d, Z1 = %d, Z2 = %d\n', L1, L2, Z1, Z2);

    % Bilinear interpolation
    psihat = wL1*wZ1*psigrid(L1,Z1) + wL1*wZ2*psigrid(L1,Z2) + ...
            wZ1*wL2*psigrid(L2,Z1) + wZ2*wL2*psigrid(L2,Z2);
    return

function [phi_c] = phi_c_monin_obukhov (x)
% --- Evaluate the Monin-Obukhov phi function for scalars at x

    if (x < 0)
        phi_c = (1 - 16 * x)^(-0.5);
    else
        phi_c = 1 + 5 * x;
    end
    return

function [phi_m] = phi_m_monin_obukhov (x)
% --- Evaluate the Monin-Obukhov phi function for momentum at x

    if (x < 0)
        phi_m = (1 - 16 * x)^(-0.25);
    else
        phi_m = 1 + 5 * x;
    end
    return

function [psi_c] = psi_c_monin_obukhov (x)
% --- Evaluate the Monin-Obukhov psi function for scalars at x

    if (x < 0)
        y = (1 - 16 * x)^0.25;
        psi_c = 2 * log((1 + y^2)/2);
    else
        psi_c = -5 * x;
    end
    return

function [psi_hat_c] = psi_c_rsl (z, h, L, c1, c2)

% --- Evaluate the roughness sublayer (RSL) function psi_hat for scalars
% at z. Note that z has already been adjusted for the displacement height
% (i.e., using z - d).

% ------------------------------------------------------
% Input
%   z            ! Vertical height - displacement height (m)
%   h            ! Canopy height - displacement height (m)
%   L            ! Obukhov length (m)
%   c1           ! Parameter for RSL function phi_hat (dimensionless)
%   c2           ! Parameter for RSL function phi_hat (dimensionless)
%
% Output
%   psi_hat_c    ! RSL psi_hat function for scalars (dimensionless)
% ------------------------------------------------------

    % The function to integrate depends on unstable (f1) or stable (f2)

    f1 = @(x) (1-16*x/L).^(-0.5) .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
    f2 = @(x) (1+5*x/L)          .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;

    % Numerically integrate the function from z to infinity

    if (L < 0)
        psi_hat_c = integral (f1, z, inf);
    else
        psi_hat_c = integral (f2, z, inf);
    end
    return

function [psi_m] = psi_m_monin_obukhov (x)
    % --- Evaluate the Monin-Obukhov psi function for momentum at x

    if (x < 0)
        y = (1 - 16 * x)^0.25;
        psi_m = 2 * log((1 + y)/2) + log((1 + y^2)/2) - 2 * atan(y) + pi / 2;
    else
        psi_m = -5 * x;
    end
    return

function [psi_hat_m] = psi_m_rsl (z, h, L, c1, c2)
% --- Evaluate the roughness sublayer (RSL) function psi_hat for momentum
% at z. Note that z has already been adjusted for the displacement height
% (i.e., using z - d).

% ------------------------------------------------------
% Input
%   z            ! Vertical height - displacement height (m)
%   h            ! Canopy height - displacement height (m)
%   L            ! Obukhov length (m)
%   c1           ! Parameter for RSL function phi_hat (dimensionless)
%   c2           ! Parameter for RSL function phi_hat (dimensionless)
%
% Output
%   psi_hat_m    ! RSL psi_hat function for momentum (dimensionless)
% ------------------------------------------------------

% The function to integrate depends on unstable (f1) or stable (f2)

    f1 = @(x) (1-16*x/L).^(-0.25) .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;
    f2 = @(x) (1+5*x/L)           .* (1 - (1 - c1*exp(-c2*x/(2*h)))) ./ x;

    % Numerically integrate the function from z to infinity

    if (L < 0)
        psi_hat_m = integral (f1, z, inf);
    else
        psi_hat_m = integral (f2, z, inf);
    end
    return







