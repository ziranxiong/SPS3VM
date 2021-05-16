clc;
clear;
global KTYPE KSCALE size_training size_data C C_star C_end s
KTYPE = 6;
KSCALE = 0.5;
s = 0;
C_star = 0.1;

C = 0.125;
C_end = 128;

norm = 'usps';
size_data = 1000;
size_training = 600;

for i=1:1
    [time,error,original_x,original_y,extended_x,extended_y,initial_solution,extended_real_y] = main(norm);
end
