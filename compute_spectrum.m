function [transmission, reflection] = compute_spectrum(wave_list, thick_list, index_refl, index_trans, index_list)

transmission = zeros(length(wave_list), 1);
reflection = transmission;

mat_in = index_refl;
mat_out = index_trans;

k_in = 2*pi*mat_in./wave_list;
k_out = 2*pi*mat_out./wave_list;

parfor i = 1:length(wave_list)
    
    wave = wave_list(i);
    
    M_i = [1, 0; 0, 1];
    k_list = 2*pi*index_list./wave;
    
    for j = 1:length(thick_list)
        kL = k_list(j)*thick_list(j);
        
        M_i = [cos(kL), 1/k_list(j)*sin(kL); -k_list(j)*sin(kL), cos(kL)]*M_i;
    end
    
    t = 2*1i*k_in(i)/(-M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2) + 1i*(k_out(i)*M_i(1,1) + k_in(i)*M_i(2,2)));
    r = ((M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2))+1i*(k_in(i)*M_i(2,2)-k_out(i)*M_i(1,1)))/ ...
        ((-M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2))+1i*(k_in(i)*M_i(2,2)+k_out(i)*M_i(1,1)))
    transmission(i) = abs(t)^2;
    reflection(i) = abs(r)^2;

end