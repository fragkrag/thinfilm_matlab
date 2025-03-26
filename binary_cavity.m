function [index_layers, thick_layers] = binary_cavity(top_dbr_num, bot_dbr_num, index_high, index_low, index_cavity, operating_wave, cavity_length)

dbr_top_index = binary_index_layers("HLH", index_high, index_low, top_dbr_num);
dbr_bot_index = binary_index_layers("HLH", index_high, index_low, bot_dbr_num);

dbr_top_thick = operating_wave./dbr_top_index/4;
dbr_bot_thick = operating_wave./dbr_bot_index/4;

index_layers = [dbr_top_index(:); index_cavity; dbr_bot_index(:)];
thick_layers = [dbr_top_thick(:); cavity_length; dbr_bot_thick(:)];

