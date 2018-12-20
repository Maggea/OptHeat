function [ SLM, rebuild ] = GerchbergSaxton( target, SLM, pup, Nit)
%GerchbergSaxton 

for it = 0 : Nit
    
    rebuild = ifftshift(fft2(fftshift(SLM)));
    
    rebuild = target .* exp( 1j .* (angle(rebuild)));
    
    SLM = pup .* exp( 1j .* angle(ifftshift(ifft2(fftshift(rebuild)))));
    
end

rebuild = ifftshift(fft2(fftshift(SLM)));  

end

