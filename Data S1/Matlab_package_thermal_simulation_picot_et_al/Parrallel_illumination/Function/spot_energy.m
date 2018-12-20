function [ spot_energy,naa,maa ] = spot_energy( slice,xcor,ycor)
%spot_energy Integrate the energy in the spot
%   Integrate the temperature in each pixel of the desired target, and then
%   divide by the amount of pixels in a comprehensive way for choosing
%   coordinates, related to figure orientation.

global dxy N spot_radius

for cran = 1 : size(xcor,2)   
    [na,ma] = meshgrid(0-N/2:N-1-N/2 ,0-N/2:N-1-N/2);    % Meshgrid for area
    maa = ma .* dxy; naa = na .* dxy;
    target = zeros(N,N);
    target = target + (hypot(naa - xcor(cran) ,maa + ycor(cran) )<= spot_radius);
    [coorX,coorY] = find(target >= 1);
    NbPixels = size(coorX,1);
    spot_energy = 0 ;
    for cranPixels = 1 : NbPixels
        spot_energy = spot_energy + slice(coorX(cranPixels),coorY(cranPixels));
    end
end

end
