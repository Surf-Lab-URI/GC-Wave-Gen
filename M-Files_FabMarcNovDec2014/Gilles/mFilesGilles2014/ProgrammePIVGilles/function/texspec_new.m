clear
data=imread('smalljump8_2000_a.bmp');
nvertstart=297;
nvertend=600;
nwiny=16;
nwinx=16;

surfdata=data(nvertstart:nvertend,:);
surfdata=double(surfdata);
% 
% [ix iy]=gradient(surfdata);
% mag=sqrt(ix.^2+iy.^2);
% template=[1 0 1
%     0 1 0
%     1 0 1];
% corrln=normxcorr2(template,mag);
% [mmm nnn]=size(corrln);
% weight=abs(corrln(2:mmm-1, 2:nnn-1).*surfdata);
% weight=uint8(weight);
% surfdata=weight;

% quantize DOWN from 8-bit (=nmax) grayscale to nquant-bit scale 
nmax=8;
nquant=4;
 
% quant=uint8(double(crap)-mod(double(crap),2^(nmax-nquant)*ones(my,mx)));
% remember that because of "ceil" the range goes from 1->2^nquant and NOT from 0->(2^nquant-1).

quant=uint8(ceil(double(surfdata)./2^(nmax-nquant))); 
quantdoub=double(quant);

[p1 p2]=size(quantdoub);

% for 0 overlap
 nwindowx=floor(p2/nwinx);
 nwindowy=floor(p1/nwiny); 
 
% for 50  overlap
%  nwindowx=2*floor(p2/nwinx)-1;
%  nwindowy=2*floor(p1/nwiny)-1;

 % for 75 overlap
 nwindowx=4*floor(p2/nwinx)-3;
 nwindowy=4*floor(p1/nwiny)-3;
 
nx=3;
ny=3; 

% for 75  overlap
for pqr=1:nwindowx;
    pqr
    vertsec=quantdoub(:,(pqr-1)*nwinx+1-(pqr-1)*3*nwinx/4:pqr*nwinx-(pqr-1)*3*nwinx/4); %50 % overlap
    
    for mkl=1:nwindowy;
        crap=vertsec((mkl-1)*nwiny+1-(mkl-1)*3*nwiny/4:mkl*nwiny-(mkl-1)*3*nwiny/4,:);  % 50 % overlap
        [my mx]=size(crap);
        
% for 50 overlap
% 
% for pqr=1:nwindowx;
%     pqr
%     vertsec=quantdoub(:,(pqr-1)*nwinx+1-(pqr-1)*nwinx/2:pqr*nwinx-(pqr-1)*nwinx/2); %50 % overlap
%     
%     for mkl=1:nwindowy;
%         crap=vertsec((mkl-1)*nwiny+1-(mkl-1)*nwiny/2:mkl*nwiny-(mkl-1)*nwiny/2,:);  % 50 % overlap
%         [my mx]=size(crap);

% for 0 overlap

% for pqr=1:nwindowx;
%     pqr
%     vertsec=quantdoub(:,(pqr-1)*nwinx+1:pqr*nwinx);
%     
%     for mkl=1:nwindowy;
%         crap=vertsec((mkl-1)*nwiny+1:mkl*nwiny,:);
%         [my mx]=size(crap);

        nxtot=floor(mx/nx);
        nytot=floor(my/ny);
        tsmatrix=0*ones(8,6561);
        hmmatrix=0*ones(8,6561);
        vmmatrix=0*ones(8,6561);
        dm1matrix=0*ones(8,6561);
        dm2matrix=0*ones(8,6561);
        kmatrix=0*ones(8,6561);
        junk=0;
        nbrhd=0;
        for j=2:my-1;
            for i=2:mx-1;
                nbrhd=crap(j-1:j+1,i-1:i+1);
                % collapsing matrix to a vector, in clockwise direction
                junk(1)=nbrhd(1,1);
                junk(2)=nbrhd(1,2);
                junk(3)=nbrhd(1,3);
                junk(4)=nbrhd(2,3);
                junk(5)=nbrhd(3,3);
                junk(6)=nbrhd(3,2);
                junk(7)=nbrhd(3,1);
                junk(8)=nbrhd(2,1);
                midjunk=nbrhd(2,2); % middle value in matrix
                % creating the texture unit
                for kk=1:nx*ny-1;
                    if(junk(kk)<midjunk)
                        ts(kk)=0;
                    elseif(junk(kk)==midjunk)
                        ts(kk)=1;
                    else
                        ts(kk)=2;
                    end
                end
                
                % ordering a = 1
                tsnew=ts;
                % horizontal structure
                if(tsnew(1)==tsnew(2)==tsnew(3))
                    pabch=3;
                elseif(tsnew(1)==tsnew(2)~=tsnew(3))|(tsnew(1)==tsnew(3)~=tsnew(2))|(tsnew(2)==tsnew(3)~=tsnew(1))
                    pabch=2;
                else
                    pabch=1;
                end
                if(tsnew(6)==tsnew(7)==tsnew(5))
                    pfghh=3;
                elseif(tsnew(6)==tsnew(7)~=tsnew(5))|(tsnew(6)==tsnew(5)~=tsnew(7))|(tsnew(7)==tsnew(5)~=tsnew(6))
                    pfghh=2;
                else
                    pfghh=1;
                end
                % vertical structure
                if(tsnew(3)==tsnew(4)==tsnew(5))
                    pabcv=3;
                elseif(tsnew(3)==tsnew(4)~=tsnew(5))|(tsnew(3)==tsnew(5)~=tsnew(4))|(tsnew(4)==tsnew(5)~=tsnew(3))
                    pabcv=2;
                else
                    pabcv=1;
                end
                if(tsnew(7)==tsnew(1)==tsnew(8))
                    pfghv=3;
                elseif(tsnew(7)==tsnew(1)~=tsnew(8))|(tsnew(7)==tsnew(8)~=tsnew(1))|(tsnew(1)==tsnew(8)~=tsnew(7))
                    pfghv=2;
                else
                    pfghv=1;
                end
                % diagonal 1 structure
                if(tsnew(1)==tsnew(4)==tsnew(2))
                    pabcd1=3;
                elseif(tsnew(1)==tsnew(4)~=tsnew(2))|(tsnew(1)==tsnew(2)~=tsnew(4))|(tsnew(2)==tsnew(4)~=tsnew(1))
                    pabcd1=2;
                else
                    pabcd1=1;
                end
                if(tsnew(7)==tsnew(5)==tsnew(8))
                    pfghd1=3;
                elseif(tsnew(7)==tsnew(5)~=tsnew(8))|(tsnew(7)==tsnew(8)~=tsnew(5))|(tsnew(5)==tsnew(8)~=tsnew(7))
                    pfghd1=2;
                else
                    pfghd1=1;
                end
                
                % diagonal 2 structure
                if(tsnew(2)==tsnew(3)==tsnew(5))
                    pabcd2=3;
                elseif(tsnew(2)==tsnew(3)~=tsnew(5))|(tsnew(2)==tsnew(5)~=tsnew(3))|(tsnew(3)==tsnew(5)~=tsnew(2))
                    pabcd2=2;
                else
                    pabcd2=1;
                end
                if(tsnew(4)==tsnew(6)==tsnew(7))
                    pfghd2=3;
                elseif(tsnew(4)==tsnew(6)~=tsnew(7))|(tsnew(4)==tsnew(7)~=tsnew(6))|(tsnew(6)==tsnew(7)~=tsnew(4))
                    pfghd2=2;
                else
                    pfghd2=1;
                end
                
                ntu=0;
                for ii=1:nx*ny-1;
                    ntu=ntu+tsnew(ii)*3^(ii-1);
                end
                tsmatrix(1,ntu+1)=tsmatrix(1,ntu+1)+1;
                hmmatrix(1,ntu+1)=pabch*pfghh;
                vmmatrix(1,ntu+1)=pabcv*pfghv;
                dm1matrix(1,ntu+1)=pabcd1*pfghd1;
                dm2matrix(1,ntu+1)=pabcd2*pfghd2;
                % central symmetry
                if (tsnew(1)==tsnew(8))&(tsnew(2)==tsnew(7))&(tsnew(3)==tsnew(6))&(tsnew(4)==tsnew(5))
                    kmatrix(1,ntu+1)=4;
                else
                    kmatrix(1,ntu+1)=0;
                end
                
                % ordering 2 - 8
                for iorder=2:8;
                    tsnew=[ts(iorder:8),tsnew(1:iorder-1)];
                    ntu=0;
                    for ii=1:nx*ny-1;
                        ntu=ntu+tsnew(ii)*3^(ii-1);
                    end
                    tsmatrix(iorder,ntu+1)=tsmatrix(iorder,ntu+1)+1;
                end
            end
        end
        tsmatrix=tsmatrix/sum(sum(tsmatrix));
        
        % black-white symmetry (independent of ordering) choosing order=a=1
        
        numer=0;
        for i=1:3280;
            numer=numer+abs(tsmatrix(1,i)-tsmatrix(1,6561-i));
        end
        bws(mkl,pqr)=(1-numer/sum(tsmatrix(1,:)))*100;
        
        % geometric symmetry
        sumgs=0.0;
        for jj=1:4
            summa=0;
            for i=1:6561
                summa=summa+abs(tsmatrix(jj,i)-tsmatrix(jj+4,i));
            end
            sumgs=sumgs+summa/(2*sum(tsmatrix(jj,:)));
        end
        gs(mkl,pqr)=(1-sumgs/4)*100;
        
        % degree of direction
        
        summm=0.0;
        for mm=1:3;
            sumnn=0.0;
            for nn=mm+1:4;
                sumii=0.0;
                for ii=1:6561;
                    sumii=sumii+abs(tsmatrix(mm,ii)-tsmatrix(nn,ii));
                end
                sumnn=sumnn+sumii/(2*sum(tsmatrix(mm,:)));
            end
            summm=summm+sumnn;
        end
        dd(mkl,pqr)=(1-summm/6)*100;
            
        % orientational features
        
        % horizontal structure
        mhstemp=0.0;
        for mm=1:6561;
            mhstemp=mhstemp+tsmatrix(1,mm)*hmmatrix(1,mm);
        end
        mhs(mkl,pqr)=mhstemp;

        % vertical structure
        mvstemp=0.0;
        for mm=1:6561;
            mvstemp=mvstemp+tsmatrix(1,mm)*vmmatrix(1,mm);
        end
        mvs(mkl,pqr)=mvstemp;
        
        % diagonal 1 structure
        mds1temp=0.0;
        for mm=1:6561;
            mds1temp=mds1temp+tsmatrix(1,mm)*dm1matrix(1,mm);
        end
        mds1(mkl,pqr)=mds1temp;
        
        % diagonal 2 structure
        mds2temp=0.0;
        for mm=1:6561;
            mds2temp=mds2temp+tsmatrix(1,mm)*dm2matrix(1,mm);
        end
        mds2(mkl,pqr)=mds2temp;
        
        % central symmetry
        cstemp=0.0;
        for mm=1:6561;
            cstemp=cstemp+tsmatrix(1,mm)*kmatrix(1,mm);
        end
        cs(mkl,pqr)=cstemp;
        
% for 75 overlap

        yaxis(mkl,pqr)=nvertstart+(mkl+1)*nwiny/4;
        xaxis(mkl,pqr)=(pqr+1)*nwinx/4;
    end
end

% 50  overlap
%         
%         yaxis(mkl,pqr)=nvertstart+(mkl-1)*nwiny+1-mkl*nwiny/2;
%         xaxis(mkl,pqr)=(pqr-1)*nwinx+1-pqr*nwinx/2;
%         
% for 0  overlap
%         
%         yaxis(mkl,pqr)=nvertstart+(mkl-1)*nwiny+nwiny/2;
%         xaxis(mkl,pqr)=(pqr-1)*nwinx+nwinx/2;

    end
end
surfhand