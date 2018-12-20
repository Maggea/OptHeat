function [ scat_array] = scatterers_array_generation()
%scatterers_matrix_generation Generates the 3D volume filled with
%scatterers
%   Randomly fills the volume defined before with a specific concentration
%   of scatterers, with a predefined size and indice of refraction.
    global N Nspheres_per_layer Nz
    
%%  Array generation

    scat_array                                   =  ones(N,N,Nz);           % Phase Mask from scatterers
parfor ijkl = 1 : Nz   
    dn=0.1;ww=2;
    xxx=cumsum(poisspdf([0:10000],Nspheres_per_layer(1)));
    NN=rand(1);
    [stam,Nspheres]=min(abs(xxx-NN));
    if(Nspheres_per_layer*Nspheres>0)
        for i=1:Nspheres;
            ix=floor(rand(1)*(N-1e-3))+1;
            iy=floor(rand(1)*(N-1e-3))+1;
            scat_array(:,:,ijkl)=scat_array(:,:,ijkl)+dn*exp(-(([1:N]-ix)/(ww^2)).^2)'*exp(-(([1:N]-iy)/(ww^2)).^2);
        end
    end
end

end

