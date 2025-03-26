function [index_layers, thick_layers] = edgepass(type, index_high, index_low, num_layers, operating_wave)

if type == "short"
    index_layers = binary_index_layers("LHL", index_high, index_low, num_layers);
    
elseif type == "long"
    index_layers = binary_index_layers("HLH", index_high, index_low, num_layers);
    
end

thick_layers = operating_wave./index_layers/4;
thick_layers(1) = thick_layers(1)/2;
thick_layers(end) = thick_layers(end)/2;
