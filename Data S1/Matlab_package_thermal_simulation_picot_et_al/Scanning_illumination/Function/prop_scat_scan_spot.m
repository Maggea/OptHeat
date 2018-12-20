function [ presource_calib,xx,yy ] = prop_scat_scan_spot( focus_depth,scat_array)

%prop_scat_scan_spot Propagation of the light field through the scatterers
%array for a defined spot target in X, Y and Z.
%   First, we define constants and generate arrays. Then we run the
%   GerchbergSaxton algorithm to retrieve the rebuild expected at the focal
%   plane for the spot required. After, we back propagate this rebuild to
%   the surface of the sample. Finally, we propagate the light field
%   through the scatterers array.
    global N dxy wvl NA dz Nz spot_radius

%% define units
    dp                                 =  2*pi/dxy/N;                      % spatial frequency pixel
    k                                  =  2*pi/wvl;                        % wavevector [inverse microns]
    NA_eff                             =  NA./1.3;                         % effective NA of lens (real NA divided by sample refractive index, 1.5 for PMMA)
    pmax                               =  2*pi/wvl*NA_eff;                 % maximal spatial frequency afforded by lens

    %% define x and p spaces
    x                                  =  linspace (-dxy*N/2,dxy*(N-1)/2,N);
    px                                 =  linspace (-dp*N/2,dp*(N-1)/2,N);
    y                                  =  x; py  =  px;
    pp                                 =  sqrt (px'.^2*(py*0+1)+(px'*0+1)*(py.^2));
    yy                                 =  (x')*(y*0+1); xx  = (x'*0+1)*(y);
    Aim                                =  zeros (N,N,size(focus_depth,2));
    zinit=focus_depth.*-1;
    amp_field = zeros(N,N,Nz);

%% Definition of the spot position in the focal plane with GS algorithm
pup = (hypot(xx,yy)/(N.*dxy/2.0)<1.0);     % Pupil definition
target = zeros(N,N);
target = target + (hypot(xx,yy ) <= spot_radius);                % WARNiNG : - X  ; + Y
SLM = ifftshift(ifft2(fftshift((target .* exp( i .* 2 .* pi .* rand(N,N))))));
SLM = pup .* exp( i .* angle(SLM));                                        %Phase distribution at SLM
[SLM,rebuild] = GerchbergSaxton( target, SLM, pup, 20);            %GS for rebuild

%% Back Propagation of light from focal plane to the surface of the volume
pA0=ifftshift (fft2 (fftshift (rebuild)));
kmap=sqrt (1-(pp/k).^2).*(pp<k);
zphase=exp (i*k*kmap*zinit); % phase factor for z propagation
pA=pA0 .* zphase;
pA=pA .* (pp<pmax);
Aout=ifftshift (ifft2 (fftshift (pA)));
working_step = Aout;

%% Propagation of light field through scatterers
zphase=exp (i*k*sqrt (1-(pp/k).^2)*dz);

for iz=1:Nz
    mask=(scat_array(:,:,iz)-1)*dz*k;
    pA=ifftshift (fft2 (fftshift (working_step.*exp(i*mask))));
    pA=pA.*zphase;
    Aout=ifftshift (ifft2 (fftshift (pA)));
    working_step = Aout;
    amp_field(:,:,iz) = Aout;                                              % Array which contains each XYZ field amplitude for each simulation
end

presource_calib = (abs (amp_field).^2);

end
