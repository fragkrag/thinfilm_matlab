clc;
close all;
clear all;
 
index_high = 2.4;
index_low = 1.45;

wave_list = (300:1500)*1e-9;

%edge pass high pass (hp) short pass
hp_layers = 12;
hp_operating_wave = 950e-9;

%edge pass low pass (lp) long pass
lp_layers = 12;
lp_operating_wave = 520e-9;

%cavity (cav)
cav_layers_top = 12;
cav_layers_bot = 12;
cav_operating_wave = 700e-9;
cav_index = index_low;
cav_thick = cav_operating_wave/2/cav_index;

%hp-lp coupling layer
hp_lp_couple_thick = 500e-9;
hp_lp_couple_index = index_low;

%lp-cav coupling layer
lp_cav_couple_thick = 500e-9;
lp_cav_couple_index = index_low;

%construct edgepass, cavities, and coupling layers
[hp_index, hp_thick] = edgepass("short", index_high, index_low, hp_layers, hp_operating_wave);
[lp_index, lp_thick] = edgepass("long", index_high, index_low, lp_layers, lp_operating_wave);
[cav_index, cav_thick] = binary_cavity(cav_layers_top, cav_layers_bot, index_high, index_low, cav_index, cav_operating_wave, cav_thick);

%optimizations
%FOM is sum over transmission over wavelength range
fom_wave_range = (wave_list > 480e-9) & (wave_list < 1000e-9) - ((wave_list > 650e-9) & (wave_list < 800e-9));
step_size = 2.5e-9;
total_iter = 300;
%optimized layers are the two coupling layers
opt1_index = length(hp_index)+1;
opt2_index = opt1_index+length(lp_index)+1;

full_index = [hp_index(:); hp_lp_couple_index; lp_index(:); lp_cav_couple_index; cav_index(:)];
full_thick = [hp_thick(:); hp_lp_couple_thick; lp_thick(:); lp_cav_couple_thick; cav_thick(:)];

hp_lp_opt = hp_lp_couple_thick;
lp_cav_opt = lp_cav_couple_thick;
opt_thick = full_thick;

for i = 1:total_iter
      
    %compute system
    [t, r, dt, dr] = compute_spectrum_with_partialderivative(wave_list, opt_thick, 1.45, 1.45, full_index, [opt1_index,opt2_index]);
    
    %compute fom
    T_spec(i, :) = abs(t).^2;
    T_fom(i) = sum(T_spec(i,:).*fom_wave_range, 2);
    
    %compute gradient
    dT = -2*real(conj(t).*dt).*fom_wave_range';
    grad = step_size*sum(dT./sqrt(sum(dT.^2,1)),1);
    
    %update
    hp_lp_opt = hp_lp_opt - grad(1);
    lp_cav_opt = lp_cav_opt - grad(2);
    opt_thick = [hp_thick(:); hp_lp_opt; lp_thick(:); lp_cav_opt; cav_thick(:)];
    
    hp_lp_list(i) = hp_lp_opt;
    lp_cav_list(i) = lp_cav_opt;
    
end

[val, loc] = min(T_fom);

figure
plot(wave_list, T_spec(1,:))
hold on
plot(wave_list, T_spec(loc,:))


figure
plot(hp_lp_list)
hold on
plot(lp_cav_list)
legend("hp-lp","lp-cav");


figure
plot(T_fom)
title("fom")

save("opt_results_4.mat","fom_wave_range","hp_lp_list","lp_cav_list","T_fom","wave_list","step_size","lp_layers","hp_layers","hp_operating_wave","lp_operating_wave","cav_layers_top","cav_layers_bot","cav_operating_wave","cav_thick");
