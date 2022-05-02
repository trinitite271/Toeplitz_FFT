close all
clear 
ns=1024*1;
A=(1:ns);
% A=[3,2,1];
% A=[1,2];
B=A+3;
C=B+3;
D=C+3;
E=D+3;
F=E+3;
G=F+3;
H=G+3;
A1=toeplitz(A);
B1=toeplitz(B);
C1=toeplitz(C);
D1=toeplitz(D);
E1=toeplitz(E);
F1=toeplitz(F);
G1=toeplitz(G);
H1=toeplitz(H);

Gz=[A1,B1,C1,D1;E1,F1,G1,H1];
Gzz=Gz+3;
Gzzz=Gzz+3;
Gz=[Gz,Gzz,Gzzz];
Gz1=Gz';
BB=Gz(1,:);
AA=Gz*BB';
K=Gz'*(AA);
% AAA=AA';
%%
pj=BB';
AA=(Gz*pj);
% tic
% Hj1=(Gz'*AA);
% toc
%%
Gz1=Gz';
nt=2;nr=2;
AAA=reshape(AA,[ns,nt]);
AAA=AAA';
% Z=zeros(1,ns*nr*nt*3);
Z=zeros(nr*nt*3,ns);
p_zero=zeros(1,ns);
tic
for zz=1:nt*nr*3
    gg_ext1=zeros(1,ns);
    for ii=1:nt
        D1=Gz1((zz-1)*ns+1,:);
        d1=D1((ii-1)*ns+1:ii*ns);
        
%         d2=AAA((ii-1)*ns+1:ii*ns);
        d2=AAA(ii,:);
        t_ext1=[d1 0 flip(d1(2:ns))];
%         p_zero=zeros(1,ns);
        p_ext1=[d2 p_zero];
        gg_ext=ifft(fft(t_ext1).*fft(p_ext1));
        gg_ext=gg_ext(1:ns);

%         Z(:,(zz-1)*ns+1:zz*ns)=Z(:,(zz-1)*ns+1:zz*ns)+gg_ext;
         Z(zz,:)=Z(zz,:)+gg_ext;
        
    end
%     Z((zz-1)*ns+1:zz*ns,:)=gg_ext1;
end
toc
tic
sixsixsix = matfft([nt,nr,ns],Gz1,AAA);
toc
% Z=W*Z';
Z=reshape(Z,nr*nt*ns*3,1);
% Q=Z-Hj1;
