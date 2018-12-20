global dxy dz N wvl NA spot_radius Nspheres_per_layer Nz
addpath('./Function/');
%% WARNING
% This simulation is highly computer intensive, and requires a high
% quantity of RAM (roughly 100 GB in this configuration at peak) to run with ease
% in the infinite media hypothesis.
% As a consequence, if you want to test it quickly and check that you own
% the necessary MATLAB toolbox, you should decrease N, Nslice, thickness 
% and Ntime to much lower values. Please note that for some for loops, 
% the parrallel toolbox is used to compute faster.

%% LASER SOURCE PARAMETERS
% Here we declare the parameters which are the most important for the 
% calculations. All the units are homogeneous to [J],[ms] and [um].
P = 11.e-6;                    % POWER AT THE OBJECTIVE FOR ONE SPOT[J/ms] --> 1mW is written 1.e-6     
spot_radius = 6;               % SPOT RADIUS [um] 
tau = 3;                       % ILLUMINATION LENGTH [ms] 
thickness = 450;               % THICKNESS OF SCATTERING TISSUE SIMULATED [um]
wvl = 1.030;                   % WAVELENGTH [um]

%% SPOT COORDINATES
% 0 is the center of the field for x and y ; 0 is surface for z.
% You need to set the coordinates for your spots in 3 arrays this way :
xcor = [-25 ; 0 ; 78 ];
ycor = [22 ; 0 ; 28 ];
zcor = [70 ; 150 ; 200 ];
depth_of_interest = zcor - thickness ./2 ;

%% THERMAL PARAMETERS
% Here we declare the thermal properties of the sample you want to test.
% These actual values are standard for gray brain matter, based on
% Yaroslavsky et al.,2002 & Yizhar et al. 2011.
ro = 1.e-15;                   % DENSITY [kg/um**3] 
C = 3650;                      % MASS HEAT CAPACITY [J/(kg*K)]                    
D = 140.0;                     % DIFFUSION COEFFICIENT [um**2/ms]
alpha = 0.00006;               % ABSORPTION COEFFICIENT [um-1]

%% OPTICAL PARAMETERS
% Here we declare the optical parameters of the microscope. 
NA = 0.8;                      % NUMERICAL APERTURE

%% DISCRETISATION
% Here we prepare the variables to describe our discretized volume.
% N, Nslice and Ntime are the amount of points
% you will have for X, Y, Z and time. If you just want to simulate a
% specific time frame, keep Ntime = 2 and set t = [0 time_of_interest_in_ms]; 
dxy = wvl./(2.*NA);            % XY RESOLUTION [um]
dz = 1;                        % Z RESOLUTION [um]
N = 1024;                      % AMOUNT of STEPS for X and Y
Nslice = 1024;                 % AMOUNT of STEPS for Z for whole volume
Nz = round(thickness ./ dz);  % AMOUNT of STEPS for Z for scattering volume
dt = 0.04;                     % TIME RESOLUTION [ms]
Ntime = 128;                   % AMOUNT of STEPS for TIME
t = (0 : (Ntime-1)) .* dt;     % TIME VECTOR ( starting at 0 )

%% MEMORY ALLOCATION & VOLUME DEFINITION
[u,v,w]=meshgrid(0-N/2:N-1-N/2 ,0-N/2:N-1-N/2,0-Nslice/2:Nslice-1-Nslice/2);
k = sqrt((u/(N*dxy)).^2 + (v/(N*dxy)).^2 + (w/(Nslice*dz)).^2);
clear u v w
extended_source = zeros(N,N,Nslice); % Blank Light Source array
presource = zeros(N,N,Nz); % Blank Light Source array
XYslice = zeros(N,N); % Blank XY array for spot temperature

%% SCATTERERS PHASE MASK
% Here we build the matrix filled with scatterers which will be used to
% simulate the scattering properties of gray brain matter. We highly
% recommend not to change the values, as they have been calibrated and
% verified experimentally in Papagiakoumou et al., 2013. The scatterers 
% will be randomly distributed in the volume, with a specific density of
% representation for a specific volume. 
Nspheres_per_layer =  dz.*N.*dxy.*N.*dxy./1000;   % Amount of scatterers per layer (0 for no scattering simulations) 
scat_array = scatterers_array_generation; display(['Scatterers Array : Done !'])

%% GET LIGHT REFERENCE
% The goal in this simulation is to provide to each spot, wherever they
% are in the volume, the same amount of energy to theorically induce the
% same action potential probability. Everything is based on the
% experimental in-vivo parameters described in the article. So with the
% next line, we will calibrate how much energy we should have in each spot,
% based on the propagation of a spot through 150 um of tissue.
presource_calib = prop_scat_holo_spot(150,scat_array,0,0);
int_intensity_spot_calib = spot_energy( presource_calib(:,:,150),0,0);
clear presource_calib
display(['Calibration propagation : Done !'])

%% LIGHT PROPAGATION
% Here we propagate above and below focal plane the rebuild through the
% scatterers to obtain the full volume describing the light propagation
% throught in-vivo or in-vitro sample. After each propagation, light power will be 
% compensated to match the reference value calculated above. This part is once
% again based on Papagiakoumou et al., 2013.
for n = 1 : size(xcor,1)
    display(['Spot propagation : ' num2str(n) ' / ' num2str(size(xcor,1)) ])
    if n == 1    
        presource = prop_scat_holo_spot(zcor(n),scat_array,xcor(n),ycor(n));
        int_intensity_spot = spot_energy( presource(:,:,zcor(n)),xcor(n),ycor(n));
        coef_compensation(n) =  (int_intensity_spot_calib ./ int_intensity_spot); 
        presource = presource .* coef_compensation(n);        
    else
        presource_n = prop_scat_holo_spot(zcor(n),scat_array,xcor(n),ycor(n));
        int_intensity_spot = spot_energy( presource_n(:,:,zcor(n)),xcor(n),ycor(n));      
        coef_compensation(n) =  (int_intensity_spot_calib ./ int_intensity_spot);
        presource = presource + (presource_n .* coef_compensation(n));                      
    end  
end
clear presource_n
display(['Spots propagation : Done !'])

%% SOURCE HANDLING
% Here we manipulate the 3D light volume to normalize it and scale it with
% the physical parameters of the simulation declared before. In the end, we
% get the source variable which contains the 3D light intensity map. We
% will then get S, the 3D Fourier Transform of the source, to multiply by
% the Fourier Transform of the Green's Function.
extended_source(:,:,((Nslice./2)+1)-(Nz./2):((Nslice./2)+1)+((Nz./2)-1)) = presource;
clear presource
P_total = P .* sum(coef_compensation); sum_source_p = sum(sum(extended_source,1),2);
gamma = (( P_total )./(sum_source_p((Nslice./2) ).*dxy.*dxy)); % normalization coefficient
light_source = extended_source .* gamma .* (alpha./(ro.*C)) ;
clear extended_source
S = ifftshift(fftn(fftshift(light_source)));   % 3D FFT of light_source
clear light_source
display(['Source Handling : Done !'])

%% HEAT DIFFUSION
% Here we perform the convolution of the source variable and the Green's
% function for each time frame. Currently, only the mean temperature rise
% in each spot is recorded at the end of the calculation, in the array
% "meanspot".
for time_step = (1:size(t,2))           
    display(['Heat diffusion timestep : ' num2str(time_step) ' / ' num2str(size(t,2)) ])
    real_tstep = t(time_step);
    if real_tstep < tau                % During illumination
        G = (1.0-exp(-4*D*(pi*k).^2*real_tstep))./(4.0*D*(pi*k).^2);
        [row,col] = find(k==0);
        G(row,col) = real_tstep;   
    else 
        G = (exp(-4*(pi*k).^2*D*(real_tstep-tau))-exp(-4*(pi*k).^2*D*real_tstep))...
            ./(4.0*D*(pi*k).^2);  % After the end of illumination
        [row,col] = find(k==0);
        G(row,col) = tau;  
    end
    XYZvolume = real(ifftshift(ifftn(fftshift(G.*S))));       % Convolution       
    for n = 1 : size(depth_of_interest,1)
        XYslice(:,:) = XYZvolume(:,:,(Nslice./2) + round(depth_of_interest(n)./dz));
        meanspot(time_step,n) = spot_mean(XYslice,xcor(n),ycor(n)); % Array with the mean temperature increase for each spot
    end
end
clear XYZvolume
display(['Heat Simulation : Done !'])

%% RESULTS PLOT
figure()
plot(t,mean(meanspot,2));
xlabel('Time [ms]');ylabel('Temperature rise [K]')
title('Mean temperature rise over the spots')

figure()
scatter(1:size(xcor,1),meanspot(round(tau./dt)+1,:));
xlabel('Spot ID');ylabel('Temperature rise [K]')
title('Temperature rise for each spots at the end of illumination')
axis([-inf inf 0 inf])
%%

