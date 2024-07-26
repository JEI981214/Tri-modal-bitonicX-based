close all; clear all; clc;

% [imagename1 imagepath1]=uigetfile('\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
% A=imread(strcat(imagepath1,imagename1)); 
% [imagename2 imagepath2]=uigetfile('\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
% G=imread(strcat(imagepath2,imagename2));  
% [imagename3 imagepath3]=uigetfile('\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
% B1=imread(strcat(imagepath3,imagename3));  
A=imread('MR-Gad.PNG');
G=imread('MR-T2.PNG');
B1=imread('SPECT-TI.PNG');
addpath Functions
%% Add noise
% Gaussian noise
%     sigma=40;
% if sigma>0 
%     v=sigma*sigma/(255*255);  % 
%     A=imnoise(A,'gaussian',0,v);
%     G=imnoise(G,'gaussian',0,v);
%     B1=imnoise(B1,'gaussian',0,v);
% %     figure;imshow(A,'border','tight');
% %     figure;imshow(G,'border','tight');
% %     figure;imshow(B1,'border','tight');
% end
%Poission niose
%sigma~=0;
%A=imnoise(A, 'poisson');
%G=imnoise(G, 'poisson');
%B1=imnoise(B1, 'poisson');
%% Fusion rules
sigma=0;

q=3;
w=19;
L=5;
t = 1;   
if size(A,3)>1
    A=rgb2gray(A);            
end

if size(G,3)>1
    G=rgb2gray(G);             
end

A = im2double(A);   
G = im2double(G);  
tic
%% Decomposition
AL = xbitonic2(A,q);
GL = xbitonic2(G,q);
AH=A-AL;
GH=G-GL;
%%
if  sigma~=0;
AH=xbitonic2(AH, 'additive');
GH=xbitonic2(GH,  'additive');
end 
                                 
F1 = MSMG_PAPCNN(AH,GH,t);

%% Energy-based fusion rule for energy layers
ALL = xbitonic2(AL,q);
GLL = xbitonic2(GL,q);
map1=abs(ALL>GLL);
F2=ALL.*map1+(1-map1).*GLL;
%% Fusion rule of CNP
AHH=AL-ALL;
GHH=GL-GLL;
CNP_times1 = CNP_test(abs(AHH),110);
CNP_times2 = CNP_test(abs(GHH),110);
map=(CNP_times1>=CNP_times2);
map=majority_consist_new(map,w);  
F3=map.*AHH+(1-map).*GHH;

%% Reconstruction
Fuse_img=F1+F2+F3;
M = Fuse_img;
figure,imshow(M,'border','tight');
%% Fusion of the middle fused image and the functional source image
B1= im2double(B1);
M= im2double(M);   
[hei, wid] = size(M);
%% Y channel transform
B_YUV=ConvertRGBtoYUV(B1);   
B1=B_YUV(:,:,1); 
%% Decomposition
ML = xbitonic2(M,q);
BL = xbitonic2(B1,q);
MH=M-ML;
BH=B1-BL;
if  sigma~=0;
BH=xbitonic2(BH, 'additive');
end
%% Fusion of detail layers(GEC)
V=2;U=1;r=0.01;
LGE1=(STO(MH).^V).*(clearity(MH,L).^U).*(MH.^r);         
LGE2=(STO(BH).^V).*(clearity(BH,L).^U).*(BH.^r); 
map=(LGE1>LGE2);
map=majority_consist_new(map,w);     
F1=map.*MH+~map.*BH;  
%% Energy-based fusion rule
map1=abs(ML>BL);
F2=ML.*map1+(1-map1).*BL;
%% Reconstruction
F=F1+F2;
%% YUV2RGB
F_YUV=zeros(hei,wid,3);
F_YUV(:,:,1)=F;
F_YUV(:,:,2)=B_YUV(:,:,2);
F_YUV(:,:,3)=B_YUV(:,:,3);
final_F=ConvertYUVtoRGB(F_YUV);     
toc
figure,imshow(final_F);
