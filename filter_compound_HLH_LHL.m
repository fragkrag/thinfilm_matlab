clc;
%close all;
clear all;

% index_list = [2, 1.45, 2, 1.45, 2, 1.45];
% thick_list = [200, 100, 200, 100, 200, 100]*1e-9;

index_list_1 = [1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45];

% index_list = [2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45,2,1.45];

thick_list_1 = 800e-9./index_list_1/4;
thick_list_1(1) = thick_list_1(1)/2;
thick_list_1(length(thick_list_1)) = thick_list_1(length(thick_list_1))/2;

index_list_2 = [2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45, 2, 1.45];
thick_list_2 = 650e-9./index_list_2/4;

full_index_list = [index_list_1(:); 2; index_list_2(:)];
full_thick_list = [thick_list_1(:); 800e-9/2; thick_list_2(:)];


% full_index_list = [index_list_2(:)];
% full_thick_list = [thick_list_2(:)];
% full_index_list = [index_list_1(:)];
% full_thick_list = [thick_list_1(:)];

%thick_list(14:end) = thick_list(14:end)*500e-9/550e-9;
wave_list = (400:1:900)*1e-9;

mat_in = 1;
mat_out = 1;

k_in = 2*pi*mat_in./wave_list;
k_out =2*pi*mat_out./wave_list;
spec = zeros([length(wave_list),2]);

parfor i = 1:length(wave_list)
    
    wave = wave_list(i);
    
    M_i = [1, 0; 0, 1];
    k_list = 2*pi*full_index_list./wave;
    
    for j = 1:length(full_thick_list)
        kL = k_list(j)*full_thick_list(j);
        
        M_i = [cos(kL), 1/k_list(j)*sin(kL); -k_list(j)*sin(kL), cos(kL)]*M_i;
    end
    
    t = 2*1i*k_in(i)/(-M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2) + 1i*(k_out(i)*M_i(1,1) + k_in(i)*M_i(2,2)));
    r = ((M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2))+1i*(k_in(i)*M_i(2,2)-k_out(i)*M_i(1,1)))/ ...
        ((-M_i(2,1)+k_in(i)*k_out(i)*M_i(1,2))+1i*(k_in(i)*M_i(2,2)+k_out(i)*M_i(1,1)))
    t2 = abs(t)^2;
    r2 = abs(r)^2;
    spec(i,:) = [t2, r2];
end

figure
plot(wave_list, spec)
legend('transmission','reflection')