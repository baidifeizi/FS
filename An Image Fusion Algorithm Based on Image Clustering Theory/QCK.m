clc;
clear;
close all;
% Reading visible light images
I=imread('visible light images');
num = size(I);
if numel(num)>2
   Ivis = rgb2gray(I);
else
   Ivis=I;
end
%Read infrared images
Iir=imread('Read infrared image');
% Read fused images
If1 =imread('fused image');
num = size(If1);
if numel(num)>2
   If = rgb2gray(If1);
else
   If=If1;
end
%evaluating indicator :Qabf,Labf,Labf
[Q1,L1,N1]=QLN(Ivis,Iir,If) ;

