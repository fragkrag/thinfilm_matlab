clc;
% close all;
% clear all;

index_high = 2.35;
index_low = 1.45;

contrast_lwp_LHL = sqrt(index_low)*sqrt(index_high);
contrast_hwp_LHL = contrast_lwp_LHL/sqrt(index_high)*index_low;
contrast_lwp_HLH = sqrt(index_low)*sqrt(index_high);
contrast_hwp_HLH = contrast_lwp_HLH/sqrt(index_low)*index_high;

wave_list = (300:1300)*1e-9;

%edge pass high pass (hp) short pass
hp_layers = 15;
hp_operating_wave = 940e-9;

%edge pass low pass (lp) long pass
lp_layers = 15;
lp_operating_wave = 520e-9;

%cavity (cav)
cav_layers_top = 5;
cav_layers_bot = 5;
cav_operating_wave = 700e-9;
cav_spacer_index = index_low;
cav_spacers = (500:5:900)*1e-9/2/cav_spacer_index;

cav_mid = cav_spacers(floor(length(cav_spacers)/2));

%hp-lp coupling layer
hp_lp_couple_thick = 500e-9;
hp_lp_couple_index = index_low;

%lp-cav coupling layer
lp_cav_couple_thick = 400e-9;
lp_cav_couple_index = index_low;

%construct edgepass, cavities, and coupling layers
[hp_index, hp_thick] = edgepass("short", index_high, index_low, hp_layers, hp_operating_wave);
[lp_index, lp_thick] = edgepass("long", index_high, index_low, lp_layers, lp_operating_wave);

all_trans = zeros(length(cav_spacers), length(wave_list));
all_refl = all_trans;



for i = 1:length(cav_spacers)
    [cav_index, cav_thick] = binary_cavity(cav_layers_top, cav_layers_bot, index_high, index_low, cav_spacer_index, cav_operating_wave, cav_spacers(i));
    
    full_index = [hp_index(:); hp_lp_couple_index; lp_index(:); lp_cav_couple_index; cav_index(:)];
    full_thick = [hp_thick(:); hp_lp_couple_thick; lp_thick(:); lp_cav_couple_thick; cav_thick(:)];
    
    [trans, refl] = compute_spectrum(wave_list, full_thick, 1.45, 1.45, full_index);
    all_trans(i,:) = trans;
    all_refl(i,:) = refl;
end

[cav_index, cav_thick] = binary_cavity(cav_layers_top, cav_layers_bot, index_high, index_low, cav_spacer_index, cav_operating_wave, cav_mid); 

[trans_cav, refl_cav] = compute_spectrum(wave_list, cav_thick, 1.45, 1.45, cav_index);
[trans_hp, refl_hp] = compute_spectrum(wave_list, hp_thick, 1.45, 1.45, hp_index);
[trans_lp, refl_lp] = compute_spectrum(wave_list, lp_thick, 1.45, 1.45, lp_index);

figure
subplot(4,1,1)
plot(wave_list*1e9, trans_hp)
subplot(4,1,2)
plot(wave_list*1e9, trans_cav)
subplot(4,1,3)
plot(wave_list*1e9, trans_lp)
subplot(4,1,4)
plot(wave_list*1e9, trans_hp.*trans_cav.*trans_lp)
% figure
% plot(wave_list, all_trans.')

figure
plot(wave_list, all_trans(1:10:end,:).')