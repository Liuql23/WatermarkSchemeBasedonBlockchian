function img=Anorld(watermark,T,key)

[m,n]=size(watermark);

% Arnoldº”√‹
if key==0
    N=n;
    %÷√¬“
    imgn=zeros(m,n);img=watermark;
    for t=1:T
        for ii=1:m
            for jj=1:n           
                x=mod((ii-1)+(jj-1),N)+1;
                y=mod((ii-1)+2*(jj-1),N)+1;        
                imgn(x,y)=img(ii,jj);                
            end
        end
        img=imgn;
    end
end

% ArnoldΩ‚√‹
if key==1
    N=n;
    %÷√¬“
    img=zeros(m,n);imgn=watermark;
    for t=1:T
        for ii=1:m
            for jj=1:n 
                x=mod(2*(ii-1)-(jj-1),N)+1;
                y=mod(-1*(ii-1)+(jj-1),N)+1;        
                img(x,y)=imgn(ii,jj);                
            end
        end
        imgn=img;
    end
end

end