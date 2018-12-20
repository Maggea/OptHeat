global dxy dz N wvl NA Nspheres_per_layer R Nz spot_radius
addpath('./Function/');
%% WARNING
% This simulation is highly computer intensive, and requires a high
% quantity of RAM (roughly 100 GB in this configuration at peak) to run with ease
% in the infinite media hypothesis.
% As a consequence, if you want to test it quickly and check that you own
% the necessary MATLAB toolbox, you should decrease N, Nslice, thickness 
% and Ntime to much lower values. Please note that for some for loops, 
% the parrallel toolbox is used to compute faster.

%% LASER PARAMETERS
% Here we declare the parameters which are the most important for the 
% calculations. All the units are homogeneous to [J],[ms] and [um].
P = 16.e-6;                    % POWER AT THE OBJECTIVE [J/ms] --> 1mW is written 1.e-6     
tau = 3;                       % ILLUMINATION LENGTH [ms] 
thickness = 100;               % THICKNESS OF SCATTERING TISSUE SIMULATED [um]
wvl = 1.030;                   % WAVELENGTH [um]
spot_radius = 0.4 ;            % SPOT RADIUS [um] 

%% OPTICAL PARAMETERS
% Here we declare the optical parameters of the microscope. 
NA = 0.8;                      % NUMERICAL APERTURE

%% THERMAL PARAMETERS
% Here we declare the thermal properties of the sample you want to test.
% These actual values are standard for gray brain matter, based on
% Yaroslavsky et al.,2002 & Yizhar et al. 2011.
ro = 1.e-15;                   % DENSITY [kg/um**3] 
C = 3650;                      % MASS HEAT CAPACITY [J/(kg*K)]                    
D = 140.0;                     % DIFFUSION COEFFICIENT [um**2/ms]
alpha = 0.00006;               % ABSORPTION COEFFICIENT [um-1]

%% DISCRETISATION
% Here we prepare the variables to describe our discretized volume.
% N, Nslice and Ntime are the amount of points
% you will have for X, Y, Z and time. If you just want to simulate a
% specific time frame, keep Ntime = 2 and set t = [0 time_of_interest_in_ms]; 
dxy = wvl./(2.*NA);            % XY RESOLUTION [um]
dz = 1;                        % Z RESOLUTION [um]
N = 1024;                      % AMOUNT of STEPS for X and Y
Nslice = 512;                 % AMOUNT of STEPS for Z for whole volume
Nz = round(thickness ./ dz);  % AMOUNT of STEPS for Z for scattering volume
Ntime = 1024;                  % AMOUNT of STEPS for Time

%% MEMORY ALLOCATION & VOLUME DEFINITION
KX = (-N/2:((N/2)-1))/(N*dxy);
KY = (-N/2:((N/2)-1))/(N*dxy);
KZ = (-Nslice/2:((Nslice/2)-1))/(Nslice*dz);
[kx,ky,kz] = meshgrid( KX, KY,KZ);
[u,v,w]=meshgrid(0-N/2:N-1-N/2 ,0-N/2:N-1-N/2,0-Nslice/2:Nslice-1-Nslice/2);
k = sqrt((u/(N*dxy)).^2 + (v/(N*dxy)).^2 + (w/(Nslice*dz)).^2);
XYslice = zeros(N,N,Ntime); % Blank XY array for spot temperature
XYslice_spiral = zeros(N,N,Ntime); % Blank XY array for spiral temperature
extended_source = zeros(N,N,Nslice); % Blank Light Source array

%% SPIRAL PARAMETERS
% Here we define all the key parameters for the spiral. Please note that
% here we use a spiral which will start on the edge, and end at the center
% of it.
depth = 50;                    % Focal plane depth [um]
R = 7.5 ;                      % Maximum spiral radius [um]
roll = 2 ;                     % Key  reconstruction parameter for spiral build up, do not change it.
Nsteps = 300;                  % AMOUNT of SPOTS for reconstruction
Ntours = 7;                    % AMOUNT of REVOLUTIONS
theta = linspace( 0, 2 * pi * Ntours, Nsteps );
r0 = linspace(0 , R , Nsteps);
X = fliplr((r0 .* cos(theta))); % X coordinates of the spots in the spiral
Y = fliplr((r0 .* sin(theta))); % Y coordinates of the spots in the spiral

%% TIME
% Here we define what will be our time vector, which is dependant of the
% amount of steps of our spiral.
tau_spot = tau ./ Nsteps;      % PULSE DURATION [ms]
dt = tau_spot ./ roll;         % TIME RESOLUTION [ms]
t = (0:Ntime-1) * dt;          % TIME range [ms]

%% SCATTERERS PHASE MASK
% Here we build the matrix filled with scatterers which will be used to
% simulate the scattering properties of gray brain matter. We highly
% recommend not to change the values, as they have been calibrated and
% verified experimentally in Papagiakoumou et al., 2013. The scatterers 
% will be randomly distributed in the volume, with a specific density of
% representation for a specific volume. 
Nspheres_per_layer =  (dz.*N.*dxy.*N.*dxy./1000);   % Amount of scatterers per layer (0 for no scattering simulations)       
scat_array = scatterers_array_generation; display(['Scatterers Array : Done !'])

%% LIGHT PROPAGATION
% Here we propagate above and below focal plane the rebuild through the
% scatterers to obtain the full volume describing the light propagation
% throught in-vivo sample or in-vitro.
presource = prop_scat_scan_spot(depth,scat_array);
display(['Spots propagation : Done !'])

%% SOURCE HANDLING
% Here we manipulate the 3D light volume to normalize it and scale it with
% the physical parameters of the simulation declared before. In the end, we
% get the source variable which contains the 3D light intensity map. We
% will then get S, the 3D Fourier Transform of the source, to multiply by
% the Fourier Transform of the Green's Function.
extended_source(:,:,((Nslice./2)+1)-(Nz./2):((Nslice./2)+1)+((Nz./2)-1)) = presource;
clear presource;
sum_source_p = sum(sum(extended_source,1),2);
gamma = (( P )./(sum_source_p((Nslice./2) ).*dxy.*dxy)); % normalization coefficient
light_source = extended_source .* gamma .* (alpha./(ro.*C)) ;
S = ifftshift(fftn(fftshift(light_source)));   % 3D FFT of light_source
clear light_source
display(['Source Handling : Done !'])

%% HEAT DIFFUSION FOR A SINGLE SPOT
% Here we perform the convolution of the source variable and the Green's
% function for each time frame for just a single spot.
for time_step = (1:size(t,2))            
    display(['Heat diffusion timestep : ' num2str(time_step) ' / ' num2str(size(t,2)) ])
    real_tstep = t(time_step);
    if real_tstep < tau_spot                % During illumination
        G = (1.0-exp(-4*D*(pi*k).^2*real_tstep))./(4.0*D*(pi*k).^2);
        [row,col] = find(k==0);
        G(row,col) = real_tstep;   
    else 
        G = (exp(-4*(pi*k).^2*D*(real_tstep-tau_spot))-exp(-4*(pi*k).^2*D*real_tstep))...
            ./(4.0*D*(pi*k).^2);  % After the end of illumination
        [row,col] = find(k==0);
        G(row,col) = tau_spot;  
    end
    XYZvolume = real(ifftshift(ifftn(fftshift(G.*S))));       % Convolution       
    XYslice(:,:,time_step) = XYZvolume(:,:,Nslice./2);
end
clear S G XYZvolume
display(['Heat Simulation : Done !'])

%% SPIRAL RE CONSTRUCTION
% Here we reconstruct the spiral thermal simulation by time and space
% shifting the thermal data of the single spot calculated before for each
% of the positions of the spiral
[kx,ky] = meshgrid( KX, KY);
for isteps = 1 : Nsteps    %  For each spots  
    roll_indice = isteps -1 ;  % At which roll indice are we ?
    XYslice_data = circshift(XYslice, roll.*roll_indice, 3);   %   Time shifting of the one spot simulation thermal results
    XYslice_data(:, :, 1:(roll*roll_indice)) = 0;     %   Clearing trash of time shifting heat data    
    parfor a = 1 : size(t,2)    %   For each time steps                
        XYslice_data(:,:,a) = ifftshift(fftn(fftshift(XYslice_data(:,:,a)))).* (exp(-1j.*2.*pi.*(kx.*X(isteps)+ky.*Y(isteps)))) ;  %  XY shift of the spot at t time step            
        XYslice_data(:,:,a) = real(ifftshift(ifftn(fftshift(XYslice_data(:,:,a)))));                     
        XYslice_spiral(:,:,a) = XYslice_spiral(:,:,a) + XYslice_data(:,:,a);   %  Data storage in XYslice_spiral                     
    end    
end
clear XYslice_data XYslice
display(['Spiral re construction : Done !'])


%%
for n = 1 : size(t,2)   
    meanspot(n) = spiral_mean(XYslice_spiral(:,:,n)); % Array with the mean temperature increase in the spiral
end
%%
figure()
plot(t,squeeze(XYslice_spiral(N/2 + 1,N/2 + 1,:)));hold on
plot(t,squeeze(meanspot));hold on
xlabel('Time [ms]')
ylabel('Temperature rise [K]')
title('Temperature rise for spiral configuration')
legend('Center of spiral','Mean over cell')

%%
