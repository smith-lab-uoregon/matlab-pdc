%{
This script takes a Sellmeier equations of different materials, extracts refractive indices at
given wavelengths, and calculates PDC phasematching and CW-pumped JSA functions
%}
c = physconst('LightSpeed');

% Set crystal material and polarization axes in crystallographic
% coordinate system.
% Tthere's two possibilities (consult function Sellmeier
% first): Either, 1 describes o wave, 2 the e wave, or 0, 1, 2 correspond
% to x, y, z in crystal coordinates.
crystal = "LBO";
QPM=1;
% BBO: Use x/z
% LT: Use x/z
% LBO: Use x or y pump
ax_s = [1,2]; 
ax_i = [2,1];
ax_p = [2,1];
% With angle tuning, the refractive index is interpolated by rotation from
% the first index axis to the second.

pmangle = .0; % Here we can adjust the crystal angle
pmangle = pi*pmangle/180.; % conversion to radians

lambda_p = 266; %nm, pump wavelength. This is CW
sigma_p = 0.1; % nm
lambda_s = lambda_p*2.0;
lambda_i = lambda_p*2.0;
L = 1e-3; %m
temp = 30; % degrees centigrade

% Grid size and resolution 
wl_int = 60; %nm how much to calculate around central wavelength of signal and idler
pts = 6000;

% Make signal and idler arrays
wl_s = linspace(lambda_s-wl_int/2.0,lambda_s+wl_int/2.0,pts); 
wl_i = linspace(lambda_i-wl_int/2.0,lambda_i+wl_int/2.0,pts);

% Generate 2d matrices from the 1d arrays
wl2d_s = repmat(wl_s,pts,1);
wl2d_i = repmat(transpose(wl_i),1,pts);

% alright, go to frequency:
nu_s = c*1./wl2d_s;
nu_i = c*1./wl2d_i;
nu_p = c*1./lambda_p; % Note that this is scalar (central frequency)

sigma_nu_p = (c*1./lambda_p) - (c*1./(lambda_p+sigma_p)); % Pump bandwidth in freq

% make pump function:
%pump = exp(-((nu_s+nu_i-nu_p)^2)/(2.*sigma_nu_p^2));
pump = exp(-((nu_s+nu_i-nu_p).^2.)/(2.*sigma_nu_p^2.));

% Plot that shit:
figure('Renderer', 'painters', 'Position', [0 400 1000 600])
subplot(2,3,1)
imagesc(wl_s,wl_i,pump)
colormap parula
set(gca,'YDir','normal');
title('Pump')
xlabel('\lambda_s [nm]')
ylabel('\lambda_i [nm]')

% Make a "hypothetical" pump array: every entry corresponds to a frequency
% that would give you exactly the corresponding signal and idler wavelength
wl2d_p = c*1./(nu_s+nu_i);

% Phasemtsmatch (note sin=cos interpolation between two axes)
k_p = 1e9*2.*pi*(cos(pmangle)*Sellmeier(wl2d_p,temp,ax_p(1),crystal)+sin(pmangle)*Sellmeier(wl2d_p,temp,ax_p(2),crystal))./wl2d_p;
%k_s = 1e9*2.*pi*Sellmeier(wl2d_s,temp,ax_s,crystal)./wl2d_s;
%k_i = 1e9*2.*pi*Sellmeier(wl2d_i,temp,ax_i,crystal)./wl2d_i;
k_s = 1e9*2.*pi*(cos(pmangle)*Sellmeier(wl2d_s,temp,ax_s(1),crystal)+sin(pmangle)*Sellmeier(wl2d_s,temp,ax_s(2),crystal))./wl2d_s;
k_i = 1e9*2.*pi*(cos(pmangle)*Sellmeier(wl2d_i,temp,ax_i(1),crystal)+sin(pmangle)*Sellmeier(wl2d_i,temp,ax_i(2),crystal))./wl2d_i;
% Look at center point of grid to calculate qpm for that point
middle = pts/2;
qpm = k_p(middle,middle) - k_s(middle,middle) - k_i(middle,middle);
pp_period = 2.*pi/qpm; % Theoretically required poling period
if QPM==0
    qpm=0.; % or not if turned off.
end
delta_k = k_p-k_s-k_i-qpm; % phasemismatch

% Phasematching function
phasematching = abs(sinc(delta_k.*L/2.));
% And plot
subplot(2,3,2)
imagesc(wl_s,wl_i,phasematching)
colormap parula
set(gca,'YDir','normal');
title('PM')
xlabel('\lambda_s [nm]')
ylabel('\lambda_i [nm]')

% Get full JSI:
JSI = pump.*phasematching;
figure(1)
subplot(2,3,3)
imagesc(wl_s,wl_i,JSI)
colormap parula
set(gca,'YDir','normal');
title('JSI')
xlabel('\lambda_s [nm]')
ylabel('\lambda_i [nm]')

% Plot marginal spectra and calculate bandwidths
marginal_s = trapz(JSI,1);
marginal_s = marginal_s ./ max(marginal_s);
marginal_i = trapz(JSI,2);
marginal_i = marginal_i ./ max(marginal_i);
% Fit them for FWHM extraction
fwhm = zeros(1,2);
f = fit(wl_s.',marginal_s.','gauss1') % gauss1 = single gauss function with scaling factor
coeffs = coeffvalues(f);
fwhm(1)=coeffs(3);
f = fit(wl_i.',transpose(marginal_i).','gauss1')
coeffs = coeffvalues(f);
fwhm(2)=coeffs(3);
% Plot them:
subplot(2,3,4)
plot(wl_s,marginal_s,wl_i,marginal_i)
title('Marginals')
xlabel('\lambda [nm]')
text(510,0.75,strcat('FWHM_i=',num2str(fwhm(2)),'nm'))
text(510,0.9,strcat('FWHM_s=',num2str(fwhm(1)),'nm'))

% Plot wavenumber curves:
subplot(2,3,5)
wl = linspace(250,550,500);
ax_kplot = [min([min(ax_s),min(ax_i),min(ax_p)]), max([max(ax_s),max(ax_i),max(ax_p)])]; % extracts the two involced ax indices. Will be identical for type-0
k1 = 1e-9*2.*pi*Sellmeier(wl,temp,ax_kplot(1),crystal)./wl;
k2 = 1e-9*2.*pi*Sellmeier(wl,temp,ax_kplot(2),crystal)./wl;
plot(wl,k1*1e9,wl,k2*1e9)
hold on
plot(wl2d_p(middle,middle),1e-9*k_p(middle,middle),'o')
plot(wl2d_s(middle,middle),1e-9*k_s(middle,middle),'o')
plot(wl2d_i(middle,middle),1e-9*k_i(middle,middle),'x')
hold off
title('Wavenumber')
xlabel('\lambda [nm]')
xlabel('k [nm^{-1}]')

% Helper functions:

% Here follow all the Sellmeier equations
function N = Sellmeier(wl, T, axis, crystal)
    %{
    Provides bulk Sellmeier equations for stuff.
    Congruently grown LiNbO3: The Sellmeier equations are from Edwards & Lawrence (ordinary)
    and Jundt (extraordinary)
    
    The equations are provided as functions with arguments
    being the wavelength given in [um] and the temperature given in [C]
    returning the squared refractive index.

    For Z-CUT LN (Propagation in x-direction)
    axes definition: x=0, y=2, z=2
    Crystal frame of reference! In the lab: x->z

    BBO: Bulk BBO Sellmeier, supposedly at room temperature (?), but not mentioned in the paper 
    taken from Dongxiang Zhang, Yufei Kong, Jing-yuan Zhang,
    Optics Communications; Volume 184, Issues 5–6, 15 October 2000, Pages 485-491

    BiBO: Taken from Umemura, Kentaro Miyata, Kiyoshi Kato (2007)

    LT (Lithium Tantalate): Provides bulk Sellmeier equations for congruently grown LiTaO3.
    The Sellmeier equations are from Abedin & Ito 1996
    Double check the crystal vs lab coordinate system! Not sure about this,
    but irrelevant for bulk crystals.

    KTP: Provides bulk Sellmeier equations for flux grown KTP.
    The equations are taken from Kato & Takaoka (2002)

    For KDP?

    LBO:  "Blue parametric generation from temperature-tuned LiB3O5", Frank Hanson and David Dick, 1991 
    %}

    wl = wl*1e-3; % convert from nm (written "532") to um (written "0.532")
    
    if crystal=="LN" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        oA1 = 4.9048;
        oA2 = 0.11775;
        oA3 = 0.21802;
        oA4 = 0.027153;
        oB1 = 2.2314e-8;
        oB2 = -2.9671e-8;
        oB3 = 2.1429e-8;
    
        eA1 = 5.35583; 
        eA2 = 0.100473;
        eA3 = 0.20692;
        eA4 = 100;
        eA5 = 11.34927;
        eA6 = 1.5334e-2;
        eB1 = 4.629e-7;
        eB2 = 3.862e-8;
        eB3 = -0.89e-8;
        eB4 = 2.657e-5;
    
        if axis==0 | axis==1
            oT = (T - 24.5) * (T + 570.50);
            N = oA1 + (oA2 + oB1 * oT) / (wl^2 - (oA3 + oB2 * oT)^2) + oB3 * oT - oA4 * wl^2;
        end
    
        if axis==2
            eT = (T - 24.5) * (T + 570.82);
            N = eA1 + eB1 * eT + (eA2 + eB2 * eT) / (wl^2 - (eA3 + eB3 * eT)^2) + (eA4 + eB4 * eT) / (wl^2 - eA5^2) - eA6 * wl^2;
        end
    end
    if crystal=="BBO" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if axis==0 | axis==1
            p = [2.7359, 0.01878, 0.01822, 0.01471, 0.0006081, 0.000067406];
        end
        if axis==2
            p = [2.3753, 0.01224, 0.01667, 0.01627, 0.0005716, 0.00006305];
        end
        
    	N =  p(1) + (p(2)./(wl.^2-p(3))) - p(4).*wl.^2 + p(5).*wl.^4 - p(6) * wl.^6;
     
    end
	if crystal == "LT" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        oA1 = 4.51224;
        oA2 = 0.0847522;
        oA3 = 0.19876;
        oA4 = -0.0239046;
        oB1 = -9.6649e-9;
        oB2 = 8.815e-8;
        oB3 = 4.25637e-8;
        
        eA1 = 4.52999;
        eA2 = 0.0844313;
        eA3 = 0.20344;
        eA4 = -0.0237909;
        eB1 = 1.72995e-7;
        eB2 = -4.7733e-7;
        eB3 = -8.31467e-8;
        
        if axis==0 | axis==1
            F = (T - 25.0) * (T + 25.0 + 546.0);
            N = oA1 + (oA2 + oB1 .* F) ./ (wl.^2 - (oA3 + oB2 .* F).^2) + (oB3 .* F) + (oA4 .* wl.^2);
        end
                  
        if axis==2
            F = (T - 25.0) * (T + 25.0 + 546.0);
            N = eA1 + (eA2 + eB1 .* F) ./ (wl.^2 - (eA3 + eB2 .* F).^2) + (eB3 .* F) + (eA4 .* wl.^2);
        end
    end
    if crystal == "BiBO" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if axis==0
            p = [3.07403,0.03231,0.03163,0.013376];
        end
        if axis==1
            p = [3.16940,0.03717,0.03483,0.01827];
        end
        if axis==2
            p = [3.6545,0.05112,0.03713,0.02261];
        end
        N = p(1) + p(2)/(wl^2-p(3)) - p(4)*wl^2;
    end
    if crystal == "KTP" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x1 = 3.29100;  y1 = 3.45018;  z1 = 4.59423;
        x2 = 0.04140;  y2 = 0.04341;  z2 = 0.06206;
        x3 = 0.03978;  y3 = 0.04597;  z3 = 0.04763;
        x4 = 9.35522;  y4 = 16.98825; z4 = 110.80672;
        x5 = 31.45571; y5 = 39.43799; z5 = 86.12171;
    
        if axis == 0
            N = x1 + x2 / (wl^2 - x3) + x4 / (wl^2 - x5);
        end
    
        if axis == 1
            N =  y1 + y2 / (wl^2 - y3) + y4 / (wl^2 - y5);
        end

        if axis == 2
            N =  z1 + z2 / (wl^2 - z3) + z4 / (wl^2 - z5);
        end  
    end
    if crystal == "KDP" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if axis == 0 | axis == 1 % no idea about these axes actually...
            p = [1.458524, 0.799459, 0.012741, 0.908108];
        end
        if axis == 2
            p = [1.423457, 0.708796, 0.012221, 0.225356];
        end
        N =  p(1) + p(2) / (1 - (p(3)/wl^2)) + p(4) / (1 - (5.0/wl^2)) ;
    end
    if crystal == "LBO" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if axis == 0
            p = [2.45768, -0.0098877, 0.026095, -0.013847];
        end
        if axis == 1
            p = [2.525, -0.017123, -0.0060517, -0.0087838];
        end
        if axis == 2
            p = [2.58488, -0.012737, 0.021414, -0.016293];
        end
        N =  p(1) + p(2) ./ (p(3) - wl.^2) + p(4)*wl.^2 ;
    end
end
   
    
    
    
    
    
    
    
    
    
    
    
    

