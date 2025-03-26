function [t, r, dt, dr] = compute_spectrum_with_partialderivative(wave_list, thick_list, index_refl, index_trans, index_list, deriv_layers)

t = zeros(length(wave_list), 1);
r = t;
dt = zeros(length(wave_list), length(deriv_layers));
dr = dt;

mat_in = index_refl;
mat_out = index_trans;

k_in = 2*pi*mat_in./wave_list;
k_out = 2*pi*mat_out./wave_list;

for i = 1:length(wave_list)
    
    wave = wave_list(i);
    
    M = [1, 0; 0, 1];
    k_list = 2*pi*index_list./wave;
    if length(deriv_layers) > 1
        dM = repmat(reshape(M, [1,2,2]), [length(deriv_layers), 1, 1]);
    else
        dM = M;
    end
    
    deriv_list = 1:length(deriv_layers);
    
    for j = 1:length(thick_list)
        kL = k_list(j)*thick_list(j);
        
        M_i = [cos(kL), 1/k_list(j)*sin(kL); -k_list(j)*sin(kL), cos(kL)];
        M = M_i*M;
        
        deriv_ind = find(deriv_layers == j);
        
        if deriv_ind
            
            dM_i = [-k_list(j)*sin(kL), cos(kL); -k_list(j)^2*cos(kL), -k_list(j)*sin(kL)];
            
            if length(deriv_list) > 1
                dM(deriv_ind,:,:) = dM_i*squeeze(dM(deriv_ind,:,:));
                others = deriv_list(deriv_list ~= deriv_ind);
                for k = 1:length(others)
                    dM(others(k), :, :) = M_i*squeeze(dM(others(k), :, :));
                end
                
            else
                dM = dM_i*dM;
            end
        else
            
            for k = 1:length(deriv_list)
                dM(k, :, :) = M_i*squeeze(dM(k, :, :)); 
            end
            
        end
    end
    
    t_denom = -M(2,1)+k_in(i)*k_out(i)*M(1,2) + 1i*(k_out(i)*M(1,1)+k_in(i)*M(2,2));
    r_denom = (-M(2,1)+k_in(i)*k_out(i)*M(1,2))+1i*(k_in(i)*M(2,2)+k_out(i)*M(1,1));
    r_numer = (M(2,1)+k_in(i)*k_out(i)*M(1,2))+1i*(k_in(i)*M(2,2)-k_out(i)*M(1,1));
    
    dt_denom = squeeze(-dM(:,2,1)+k_in(i)*k_out(i)*dM(:,1,2) + 1i*(k_out(i)*dM(:,1,1) + k_in(i)*dM(:,2,2)));
    
    t(i) = 2*1i*k_in(i)/(t_denom);
    r(i) = (r_numer)/(r_denom);
    dt(i,:) = conj(t(i))*2*1i*k_in(i)*dt_denom/(t_denom)^2;
    dr(i,:) = 0;
end