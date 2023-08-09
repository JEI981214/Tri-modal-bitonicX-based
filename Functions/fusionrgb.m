function F = fusionrgb(M,B1)
%FUSIONRGB 此处显示有关此函数的摘要
%% 分解
q=3;
L=5;
w=19;

ML = xbitonic2(M,q);
BL = xbitonic2(B1,q);

MH=M-ML;
BH=B1-BL;

%V=2;U=1;r=0.01;
V=2;U=1;r=0.01;
LGE1=(STO(MH).^V).*(myhighfusion2(MH,L).^U).*(MH.^r);         
LGE2=(STO(BH).^V).*(myhighfusion2(BH,L).^U).*(BH.^r); 
LGE1=(myhighfusion2(MH,L).^U).*(MH.^r);         
LGE2=(myhighfusion2(BH,L).^U).*(BH.^r); 

map=(LGE1>LGE2);
map=majority_consist_new(map,w);        
F1=map.*MH+~map.*BH;  

%% 低频融合规则
map1=abs(ML>BL);
F2=ML.*map1+(1-map1).*BL;
%% 初始融合图像
F=F1+F2;

end

