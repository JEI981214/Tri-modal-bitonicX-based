function [outputArg1] = clearity(I1,L)

I1=randn(256,256);
I=I1;
d=L;
[m n] = size(I);
%FPCM=Lo_FeaturePCM(I); %计算相位一致性
newI = ReflectEdge(I,d); %对原图进行扩展
%FPCM = ReflectEdge(FPCM,d);
resultI = zeros(m,n);
width = 2*d+1;
%Distance
% D = fspecial('gaussian',[width,width],sigma(1));
S = zeros(width,width); %pix Similarity
%h = waitbar(0,'Applying bilateral filter...');
%set(h,'Name','Bilateral Filter Progress');
for i=1+d:m+d
    for j=1+d:n+d
        pixValue =newI(i-d:i+d,j-d:j+d);
        %FPCMM=abs(FPCM(i-d:i+d,j-d:j+d));
        subValue = (abs(pixValue-newI(i,j))).^2;%求清晰度
        S =subValue;
        H = S;
%         resultI(i-d,j-d) = sum(abs(pixValue(:)).*H(:).*(FPCMM(:).^2));
          %resultI(i-d,j-d) = sum(abs(pixValue(:)).*H(:).*(FPCMM(:).^2)); 
          resultI(i-d,j-d) = sum(abs(pixValue(:)).*H(:)); 
    end
  
end
outputArg1 = resultI;
end

