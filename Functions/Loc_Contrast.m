function Y=Loc_Contrast(matrix1,matrix2)
% Compute the local contrast in NSCT domain
[row, column]=size(matrix1);
window_wide=5;
spread=(window_wide-1)/2;
cp=zeros(row,column);
matrix_en1=padarray(matrix1,[spread spread],'symmetric');
for i=1:row
    for j=1:column
        window=matrix_en1(i:1:(i+2*spread),j:1:(j+2*spread));
        cp(i,j)=sum(window(:))/(window_wide^2);
         Y(i,j)=abs(matrix2(i,j))/(cp(i,j)+0.000001);
    end
end
