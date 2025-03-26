function [t, r, dt, dr] = compute_spectrum_with_derivative(wave_list, thick_list, index_refl, index_trans, index_list, index_derivative)

t = zeros(length(wave_list), 1);
r = t;
dt = t;
dr = t;

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
        if j == index_derivative
            dM_i = [-k_list(j).*sin(kL), cos(kL); -k_list(j)^2*cos(kL), -k_list(j)*sin(kL)]*M_i;
        elseif j > index_derivative
            dM_i = [cos(kL), 1/k_list(j)*sin(kL); -k_list(j)*sin(kL), cos(kL)]*dM_i;
        end
    end
    
    t_denom = -M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2) + 1i*(k_out(i)*M_i(1,1)+k_in(i)*M_i(2,2));
    r_denom = (-M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2))+1i*(k_in(i)*M_i(2,2)+k_out(i)*M_i(1,1));
    r_numer = (M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2))+1i*(k_in(i)*M_i(2,2)-k_out(i)*M_i(1,1));
    
    dt_denom = -dM_i(2,1)+k_in(i)*k_out(i)*dM_i(1,2) + 1i*(k_out(i)*dM_i(1,1) + k_in(i)*dM_i(2,2));
    
    t(i) = 2*1i*k_in(i)/(t_denom);
    r(i) = (r_numer)/(r_denom);
    dt(i) = conj(t)*2*1i*k_in(i)*dt_denom/(t_denom)^2;
    dr(i) = 0;
end