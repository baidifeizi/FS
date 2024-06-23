clear
close all
%Read visible light images and infrared images
[filename,pathname]=uigetfile({'*.jpg;*.bmp;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'});
I1 =( imread([pathname,filename]));

[filename,pathname]=uigetfile({'*.jpg;*.bmp;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'});
I2 =( imread([pathname,filename]));

tic
Iy=histeq(I1);
V1=std2(I1);
V2=std2(Iy);
I=((V1/(V1+V2))*Iy+(V2/(V1+V2))*(I1));
%%Visible light image processing

VI1=std2(uint8(I));
VI2=std2(uint8(I2));


f_ori=I1;
[m,n]=size(f_ori);
f_ori = imresize(f_ori, [m/10, n/10], 'bilinear');

[rows,cols,dim]=size(f_ori); 
if dim==3 
F=rgb2lab(f_ori); 
else
F=double(f_ori); 
end
figure,imshow(f_ori) 
Data=reshape(F,rows*cols,dim);
cluster_n=2; 
gama=0.2; 
maxIter=50;
%%Multi scale feature extraction
[outA,outB,outObj,outNumIter]=RSSFCA(Data',gama,maxIter,cluster_n); 
[~,L]=max(outA);
Label=reshape(L,rows,cols);
Th=cluster_n; 
RL=CCF_ADB(Label,Th); 
R9=Label_image(f_ori,RL); 
R8 = imresize(R9, [m, n], 'bilinear');
%morphological heterogeneity processing 
se=strel('disk',350);
A1=imerode(R8,se); 
A1=imdilate(A1,se);

%Infrared image processing
f_ori=I2;
f_ori = imresize(f_ori, [m/10, n/10], 'bilinear');
[rows,cols,dim]=size(f_ori); 
if dim==3 
F=rgb2lab(f_ori); 
else
F=double(f_ori); 
end
figure,imshow(f_ori) 
Data=reshape(F,rows*cols,dim); 
cluster_n=2; 
gama=0.2; 
maxIter=50; 
[outA,outB,outObj,outNumIter]=RSSFCA(Data',gama,maxIter,cluster_n); 
[~,L]=max(outA);
Label=reshape(L,rows,cols);
Th=cluster_n; 
RL=CCF_ADB(Label,Th);
R9=Label_image(f_ori,RL); 
R9 = imresize(R9, [m, n], 'bilinear');
%morphological heterogeneity processing 
se=strel('disk',300);
A2=imerode(R9,se); 
A2=imdilate(A2,se);


D1 = imsubtract(double(I), double(A1)); % Bright details of visible light
L1 = imsubtract(double(A1), double(I)); % The dark details of visible light

D2 = imsubtract(double(I2), double(A2));%The bright details of infrared radiation
L2 = imsubtract(double(A2), double(I2)); %The dark details of infrared radiation

%Fusion strategy
%Processing of underlying layer information
H2=min(A1,A2);

 Aa=max((I),(I2));
   Ai=min((I),(I2));

   H2=imadd(double(H2),double(Aa)*0.1);
   H2=imsubtract(H2,double(Ai)*0.1);
 
%High frequency detail information processing
D_fused = D1+D2; 
L_fused = min(L1,L2); 

threshold = 5;
L_fused(L_fused < threshold) = 0; 

A = imadd(double(H2), D_fused); 

A = imsubtract(A, L_fused); 

A = max(min(A, 255), 0); 
A = uint8(A); 

figure('Name', 'reslut;'),imshow(A);

