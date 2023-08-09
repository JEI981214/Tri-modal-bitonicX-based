function Fire_series=CNP_test(I,N)

    theata0=0.3;
    [m,n]=size(I);
    Y=zeros(m,n);
    Y=double(Y);
    L=zeros(m,n);
    L=double(L);
    F=zeros(m,n);
    F=double(F);
    theata=theata0.*ones(m,n);
    p=zeros(m,n);
    Fire_series=zeros(m,n);
    S=I;
    T=zeros(m,n);
    link_arrange=7;
    center_x=round(link_arrange/2);
    center_y=round(link_arrange/2);
    W=zeros(link_arrange,link_arrange);
    for i=1:link_arrange
        for j=1:link_arrange
            if(i==center_x)&&(j==center_y)
                W(i,j)=0;
            else
                W(i,j)=1./sqrt((i-center_x).^2+(j-center_y).^2);
            end
        end
    end

    for t=1:N
        work=conv2(p,W,'same');
        for row=1:m
            for col=1:n
                if p(row,col)==0
                    F(row,col)=F(row,col)+S(row,col)+work(row,col);
                    L(row,col)=L(row,col)+work(row,col);
                    theata(row,col)=theata(row,col);
                else
                    F(row,col)=S(row,col)+work(row,col);
                    L(row,col)=work(row,col);
                    theata(row,col)=p(row,col);

                end
            end
        end
        U=F.*(1+L);
        p=(U>=theata);
        p=double(p);
        Fire_series=Fire_series+p;
    end

end
