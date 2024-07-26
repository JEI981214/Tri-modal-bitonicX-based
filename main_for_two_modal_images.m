close all; clear all; clc;

% [imagename1 imagepath1]=uigetfile('\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
% A=imread(strcat(imagepath1,imagename1)); 
% [imagename2 imagepath2]=uigetfile('\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
% G=imread(strcat(imagepath2,imagename2));  
% [imagename3 imagepath3]=uigetfile('\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
% B1=imread(strcat(imagepath3,imagename3));  

G=imread('MR-T2.PNG');% gray
B1=imread('SPECT-TI.PNG');%color
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

if size(G,3)>1
    G=rgb2gray(G);             
end

G = im2double(G);  
tic
%% Decomposition

%% Fusion of the middle fused image and the functional source image
B1= im2double(B1);
G= im2double(G);   
[hei, wid] = size(G);
%% Y channel transform
B_YUV=ConvertRGBtoYUV(B1);   
B1=B_YUV(:,:,1); 
%% Decomposition
GL = xbitonic2(G,q);
BL = xbitonic2(B1,q);
GH=G-GL;
BH=B1-BL;
if  sigma~=0;
BH=xbitonic2(BH, 'additive');
end
%% Fusion of detail layers(GEC)
V=2;U=1;r=0.01;
LGE1=(STO(GH).^V).*(clearity(GH,L).^U).*(GH.^r);         
LGE2=(STO(BH).^V).*(clearity(BH,L).^U).*(BH.^r); 
map=(LGE1>LGE2);
map=majority_consist_new(map,w);     
F1=map.*GH+~map.*BH;  
%% Energy-based fusion rule
map1=abs(GL>BL);
F2=GL.*map1+(1-map1).*BL;
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
