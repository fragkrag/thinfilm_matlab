clc;
close all;
clear all;

index_high = 2.4;
index_low = 1.45;

index1 = binary_index_layers("HLH", index_high, index_low, 3);

index2 = binary_index_layers("HLH", index_high, index_low, 5);
index3 = binary_index_layers("LHL", index_high, index_low, 3);

index4 = binary_index_layers("LHL", index_high, index_low, 5);