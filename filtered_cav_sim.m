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
hp_lp_couple_thick = 400e-9;
hp_lp_couple_index = index_low;

%lp-cav coupling layer
lp_cav_couple_thick = 400e-9;
lp_cav_couple_index = index_low;

%construct edgepass, cavities, and coupling layers
[hp_index, hp_thick] = edgepass("short", index_high, index_low, hp_layers, hp_operating_wave);
[lp_index, lp_thick] = edgepass("long", index_high, index_low, lp_layers, lp_operating_wave);
[cav_index, cav_thick] = binary_cavity(cav_layers_top, cav_layers_bot, index_high, index_low, cav_index, cav_operating_wave, cav_thick);

full_index = [hp_index(:); hp_lp_couple_index; lp_index(:); lp_cav_couple_index; cav_index(:)];
full_thick = [hp_thick(:); hp_lp_couple_thick; lp_thick(:); lp_cav_couple_thick; cav_thick(:)];

[trans, refl] = compute_spectrum(wave_list, full_thick, 1.45, 1, full_index);

[trans_hp, refl_hp] = compute_spectrum(wave_list, hp_thick, 1.45, 1, hp_index);
[trans_lp, refl_lp] = compute_spectrum(wave_list, lp_thick, 1.45, 1, lp_index);
[trans_cav, refl_cav] = compute_spectrum(wave_list, cav_thick, 1.45, 1, cav_index);

figure
plot(wave_list, trans)
hold
plot(wave_list, refl)
legend("trans","refl")
title('full')

figure
plot(wave_list, trans_cav+1)
hold on
plot(wave_list, trans_hp)
plot(wave_list, trans_lp+2)
legend("cav","hp","lp")
title('trans')

figure
plot(wave_list, refl_cav+1)
hold on
plot(wave_list, refl_hp)
plot(wave_list, refl_lp+2)
legend("cav","hp","lp")
title('refl')