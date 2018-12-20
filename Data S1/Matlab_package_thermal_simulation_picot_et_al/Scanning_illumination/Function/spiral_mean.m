function [ spot_mean,naa,maa ] = spiral_mean( slice)
%sumonospot_holo Compute the mean Temperature value over spots
%   Integrate the temperature in each pixel of the desired target, and then
%   divide by the amount of pixels in a comprehensive way for choosing
%   coordinates, related to figure orientation.

global dxy N R

[na,ma] = meshgrid(0-N/2:N-1-N/2 ,0-N/2:N-1-N/2);    % Meshgrid for area
maa = ma .* dxy; naa = na .* dxy;
target = zeros(N,N);
target = target + (hypot(naa ,maa  )<= R);
[coorX,coorY] = find(target >= 1);
NbPixels = size(coorX,1);
spot_mean = 0 ;
for cranPixels = 1 : NbPixels
    spot_mean = spot_mean + slice(coorX(cranPixels),coorY(cranPixels));
end
spot_mean = spot_mean ./ NbPixels;


end