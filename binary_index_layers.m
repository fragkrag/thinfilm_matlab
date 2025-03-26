function layers = binary_index_layers(stack, index_high, index_low, num_stacks)

if stack == "HLH"
    index_1 = index_high;
    index_2 = index_low;
elseif stack == "LHL"
    index_1 = index_low;
    index_2 = index_high;
end

layers = zeros(num_stacks*2+1, 1);
layers(1) = index_1;

for i = 1:(num_stacks*2)
    
    if mod(i,2) == 0
    
        layers(i+1) = index_1;
        
    else
        
        layers(i+1) = index_2;
        
    end
    
end
    