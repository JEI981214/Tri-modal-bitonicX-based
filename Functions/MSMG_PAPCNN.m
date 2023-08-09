function cp=MSMG_PAPCNN(matrixA,matrixB,t)

%%
MSMG_A=multiscale_morph(abs(matrixA),t);
MSMG_B=multiscale_morph(abs(matrixB),t);


%%
matrixA1=abs(matrixA);
matrixB1=abs(matrixB);
Para1.iterTimes=110;
Para2.iterTimes=110;
stt1=std2(matrixA1);
stt2=std2(matrixB1);

S21=max(max(matrixA1));
S22=max(max(matrixB1));

S11=graythresh(matrixA1);
S12=graythresh(matrixB1);
Para1.alpha_f=log10(1.0/stt1);
Para2.alpha_f=log10(1.0/stt2);
Para1.lambda=(S21/S11-1.0)/6;
Para2.lambda=(S22/S12-1.0)/6;
Para1.V_e=exp(-Para1.alpha_f)+S21/S11;
Para2.V_e=exp(-Para2.alpha_f)+S22/S12;
m1=(1-exp(-3* Para1.alpha_f))/(1- exp(-Para1.alpha_f))+(S21/S11-1)*exp(-Para1.alpha_f);
m2=(1-exp(-3* Para2.alpha_f))/(1- exp(-Para2.alpha_f))+(S22/S12-1)*exp(-Para2.alpha_f);
Para1.alpha_e=log( Para1.V_e / (S11 * m1));
Para2.alpha_e=log( Para2.V_e / (S12 * m2));



PCNN_timesA= PA_PCNN(MSMG_A,Para1);
PCNN_timesB= PA_PCNN(MSMG_B,Para2);
% PCNN_timesA= PA_PCNN(matrixA,Para1);
% PCNN_timesB= PA_PCNN(matrixB,Para2);

map=(PCNN_timesA>=PCNN_timesB);
%%

% matrixA=xbitonic2(matrixA,5);
% matrixB=xbitonic2(matrixB,5);

%%
cp=map.*matrixA+~map.*matrixB;