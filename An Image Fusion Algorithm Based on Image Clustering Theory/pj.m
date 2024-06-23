clc;
clear;
close all;
%%%%%%%%%%%%%%OUR

%Experimental evaluation
I1=imread('Input fused image');
H1=hx(I1);
 [im1,iv1] = variance(I1);
SF1=space_frequency(I1);
AG1=avg_gradient(I1);
EI1= edge_intensity(I1);
disp('AG:');
disp(AG1);
disp('H:');
disp(H1);
disp('SD:');
disp(iv1);
disp('SF:');
disp(SF1);

disp('EI:');
disp(EI1);