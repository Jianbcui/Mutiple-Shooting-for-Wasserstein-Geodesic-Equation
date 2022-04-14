 clear all
 
a=tic;
T=1;
% small interval numbers
%K=10;
N=20;
K=40;
tt=[0:1/K:1];
%xL=-0.5; xR=2.5;
xL=0; xR=2;
dx=1/64;
x=[xL:dx:xR];
M=length(x);
epsilon=10^(-10);
r=0.0001;

% p0(:,1)=exp(-15*(x-0.4).^2)+r;
% B=sum(p0(:,1));
% p0(:,1)=p0(:,1)/B/dx;
% p0(:,1)=ones(M,1);
% B=sum(p0(:,1));
% p0(:,1)=p0(:,1)/B/dx;
% p0(:,1)=exp(-50*(x-0.4).^2)+r;
% B=sum(p0(:,1));
% p0(:,1)=p0(:,1)/B/dx;
 %p0(:,1)=exp(-15*(x-0.4).^2)+r;
 p0(:,1)=ones(M,1);
 B=sum(p0(:,1));
 p0(:,1)=p0(:,1)/B/dx;


% lambda=0.1; for multiple shooting 


lambda=1;
 p0(:,K+1)=lambda*exp(-25*lambda*(x-1).^2)+r;
       B=sum(p0(:,K+1));
p0(:,K+1)=p0(:,K+1)/B/dx;


% 
%  p0(:,K+1)=exp(-15*(x-0.4-1*lambda).^2)+r;
%        B=sum(p0(:,K+1));
% p0(:,K+1)=p0(:,K+1)/B/dx;
% % 
 for j=1:K-1
     p0(:,j+1)=p0(:,1)+j*(p0(:,K+1)-p0(:,1))/K*lambda;
 end
 
     p0(:,K+1)=p0(:,1)+(p0(:,K+1)-p0(:,1))*lambda;

% for j=1:K-1
% %      p0(:,j+1)=exp(-15*(x-0.4).^2-lambda*(15*(x-1.5).^2-15*(x-0.4).^2)*j/K)+r;
%         p0(:,j+1)=exp(-15*(x-0.4-j/K*lambda*(1-0.4)).^2)+r;
%      B=sum(p0(:,j+1));
%      p0(:,j+1)=p0(:,j+1)/B/dx;    
% end


% for j=1:K-1
%      p0(:,j+1)=exp(-15*(x-0.4-1*lambda*j/K).^2)+r;
%      B=sum(p0(:,j+1));
%      p0(:,j+1)=p0(:,j+1)/B/dx;    
% end

for l=1:K 
V0(1:M-1,l)=ones(M-1,1)*10^-3;
end

for l=1:K+1
    pp0(1:M-1,l)=p0(1:M-1,l);
end

Y(:,1)=V0(:,1);
for i=1:K-1
  % p
  Y(:,2*i)=pp0(:,i+1);
  % v
  Y(:,2*i+1)=V0(:,i+1);
end 

 Yy(1:M-1)=Y(:,1);
for i=2:2*K-1
 Yy((M-1)*(i-1)+1:i*(M-1))=Y(:,i);
end

YY(:,1)=zeros((M-1),1);
AP=zeros((M-1),K-1);
AV=zeros((M-1),K-1);
[YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde1dmatrix(tt(1),tt(2),N,V0(:,1),pp0(:,1),M);
%mshootpde1dmatrixpercen
for i=1:K-1
    CV(:,i)=Y(:,2*i+1);
    CP(:,i)=Y(:,2*i);
end 
Tt=zeros(2,K-1);
for i=1:K-1
    Tt(:,i)=[tt(i+1),tt(i+2)];
end 
parfor i=1:K-1
    [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde1dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
end
   
for i=2:K
    YY(:,2*i)=AV(:,i-1);
    YY(:,2*i+1)=AP(:,i-1);
end 

for i=1:K-1
    %p 
F((2*(i-1))*(M-1)+1:(2*(i-1)+1)*(M-1))= YY(:,2*i+1)-Y(:,2*i);
    %v
F((2*i-1)*(M-1)+1:(2*i)*(M-1))= YY(:,2*i)-Y(:,2*i+1);
end
F((2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1)))=YY(:,2*K+1)-pp0(1:M-1,K+1);


A(1:M-1,1:M-1)=VV(1:M-1,M:2*(M-1),1); % partial p1 \partial v0
A(1:M-1,M:2*(M-1))=-eye(M-1,M-1); % partial p1 \partial v1
A(M:2*(M-1),1:M-1)=VV(M:2*(M-1),M:2*(M-1),1); %\partial v1 \partial v0
A(M:2*(M-1),2*(M-1)+1:3*(M-1))=-eye(M-1,M-1); %\partial v1 \partial p1
for i=2:K-1
    A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)-1)*(M-1)+1:(2*(i-1)+1)*(M-1))=V(:,:,i-1);
    A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)+1)*(M-1)+1:(2*(i)+1)*(M-1))=-eye(2*(M-1),2*(M-1));
end

A(2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1),(2*(K-2)+1)*(M-1)+1:(2*(K-1)+1)*(M-1))=V(1:M-1,:,K-1);
%A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
A=sparse(A);
YYY=lsqminnorm(-A,F',1.e-6);
% YYY=-A\F';
% Yy=Yy+YYY';

YYy=Yy+YYY';
Yy=YYy;
m=1
norm(F,inf)


while (norm(F,inf)/max(1,norm(Yy,inf))>10^(-5) && m<20)
 
  Y(:,1)=Yy(1:M-1);
for i=2:2*K-1
 Y(:,i)=Yy((M-1)*(i-1)+1:i*(M-1));
end

YY(:,1)=zeros((M-1),1);
AP=zeros((M-1),K-1);
AV=zeros((M-1),K-1);
[YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde1dmatrix(tt(1),tt(2),N,Y(:,1),pp0(:,1),M);

for i=1:K-1
    CV(:,i)=Y(:,2*i+1);
    CP(:,i)=Y(:,2*i);
end 
Tt=zeros(2,K-1);
for i=1:K-1
    Tt(:,i)=[tt(i+1),tt(i+2)];
end 
parfor i=1:K-1
    [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde1dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
end
   
for i=2:K
    YY(:,2*i)=AV(:,i-1);
    YY(:,2*i+1)=AP(:,i-1);
end 

for i=1:K-1
    %p 
F((2*(i-1))*(M-1)+1:(2*(i-1)+1)*(M-1))= YY(:,2*i+1)-Y(:,2*i);
    %v
F((2*i-1)*(M-1)+1:(2*i)*(M-1))= YY(:,2*i)-Y(:,2*i+1);
end
F((2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1)))=YY(:,2*K+1)-pp0(1:M-1,K+1);


A(1:M-1,1:M-1)=VV(1:M-1,M:2*(M-1),1); % partial p1 \partial v0
A(1:M-1,M:2*(M-1))=-eye(M-1,M-1); % partial p1 \partial v1
A(M:2*(M-1),1:M-1)=VV(M:2*(M-1),M:2*(M-1),1); %\partial v1 \partial v0
A(M:2*(M-1),2*(M-1)+1:3*(M-1))=-eye(M-1,M-1); %\partial v1 \partial p1
for i=2:K-1
    A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)-1)*(M-1)+1:(2*(i-1)+1)*(M-1))=V(:,:,i-1);
    A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)+1)*(M-1)+1:(2*(i)+1)*(M-1))=-eye(2*(M-1),2*(M-1));
end

A(2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1),(2*(K-2)+1)*(M-1)+1:(2*(K-1)+1)*(M-1))=V(1:M-1,:,K-1);
%A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
A=sparse(A);
YYY=lsqminnorm(-A,F',1.e-6);
% YYY=-A\F';
% Yy=Yy+YYY';

YYy=Yy+YYY';
Yy=YYy;

norm(F, inf)
m=m+1;
toc(a)
end


 Y(:,1)=Yy(1:M-1);
for i=2:2*K-1
 Y(:,i)=Yy((M-1)*(i-1)+1:i*(M-1));
end

vv0=(Y(1:M-1))';
for i=1:K-1
  % p
  pp(:,i)=Y((2*i-1)*(M-1)+1:2*i*(M-1));
  % S
  VVV(:,i)=Y(2*i*(M-1)+1:(2*i+1)*(M-1));
end   
% SS;
% pp;

for i=1:K-1
   pmed(1:M-1,i)=Y(:,2*i);
    end 
   pmed(M,i)=(1-sum(Y(:,2*i))*dx)/(dx);
 ppp(:,2:K)=pmed;  
 ppp(:,1)=p0(:,1);
 ppp(:,K+1)=p0(:,K+1);
 


% % first step finished if using multiple shooting 
%     lambda=0.1;
% % second step repeat the initial guess 
% 
% for q=1:10
%     lambda=lambda+0.9/10
%     
% V0(:,1)=(Yy(1:M-1))';
% for i=1:K-1
%   % p
%   pp0(1:M-1,i+1)=Yy((2*i-1)*(M-1)+1:2*i*(M-1));
%   % S
%   V0(:,i+1)=Yy(2*i*(M-1)+1:(2*i+1)*(M-1));
% end   
% 
% Y(:,1)=V0(:,1);
% for i=1:K-1
%   % p
%   Y(:,2*i)=pp0(1:M-1,i+1);
%   % v
%   Y(:,2*i+1)=V0(:,i+1);
% end 
% %  p0(:,K+1)=exp(-15*(x-0.4-1*lambda).^2)+r;
% %        B=sum(p0(:,K+1));
% % p0(:,K+1)=p0(:,K+1)/B/dx;
% %  p0(:,K+1)=exp(-15*(x-0.4).^2-lambda*(15*(x-1.5).^2-15*(x-0.4).^2))+r;
% %        B=sum(p0(:,K+1));
% % p0(:,K+1)=p0(:,K+1)/B/dx;
%  p0(:,K+1)=lambda*exp(-25*lambda*(x-1).^2)+r;
%        B=sum(p0(:,K+1));
% p0(:,K+1)=p0(:,K+1)/B/dx;
% 
% % p0(:,K+1)=p0(:,1)+(pin0'-p0(:,1))*lambda;
% % p0(:,K+1)=exp(-15*(x-0.4-1*lambda).^2)+r;
% %        B=sum(p0(:,K+1));
% % p0(:,K+1)=p0(:,K+1)/B/dx;
% %  p0(:,K+1)=exp(-50*(x-0.4-1*lambda).^2)+r;
% %        B=sum(p0(:,K+1));
% % p0(:,K+1)=p0(:,K+1)/B/dx;
% %  p0(:,K+1)=exp(-15*(x-1.4+1*lambda).^2)+r;
% %        B=sum(p0(:,K+1));
% % p0(:,K+1)=p0(:,K+1)/B/dx;
% 
% 
% 
% 
% 
%  Yy(1:M-1)=Y(:,1);
% for i=2:2*K-1
%  Yy((M-1)*(i-1)+1:i*(M-1))=Y(:,i);
% end
% 
% 
% 
% YY(:,1)=zeros((M-1),1);
% AP=zeros((M-1),K-1);
% AV=zeros((M-1),K-1);
% [YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde1dmatrix(tt(1),tt(2),N,V0(:,1),pp0(:,1),M);
% 
% for i=1:K-1
%     CV(:,i)=Y(:,2*i+1);
%     CP(:,i)=Y(:,2*i);
% end 
% Tt=zeros(2,K-1);
% for i=1:K-1
%     Tt(:,i)=[tt(i+1),tt(i+2)];
% end 
% parfor i=1:K-1
%     [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde1dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
% end
%    
% for i=2:K
%     YY(:,2*i)=AV(:,i-1);
%     YY(:,2*i+1)=AP(:,i-1);
% end 
% 
% for i=1:K-1
%     %p 
% F((2*(i-1))*(M-1)+1:(2*(i-1)+1)*(M-1))= YY(:,2*i+1)-Y(:,2*i);
%     %v
% F((2*i-1)*(M-1)+1:(2*i)*(M-1))= YY(:,2*i)-Y(:,2*i+1);
% end
% F((2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1)))=YY(:,2*K+1)-p0(1:M-1,K+1);
% 
% 
% A(1:M-1,1:M-1)=VV(1:M-1,M:2*(M-1),1); % partial p1 \partial v0
% A(1:M-1,M:2*(M-1))=-eye(M-1,M-1); % partial p1 \partial v1
% A(M:2*(M-1),1:M-1)=VV(M:2*(M-1),M:2*(M-1),1); %\partial v1 \partial v0
% A(M:2*(M-1),2*(M-1)+1:3*(M-1))=-eye(M-1,M-1); %\partial v1 \partial p1
% for i=2:K-1
%     A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)-1)*(M-1)+1:(2*(i-1)+1)*(M-1))=V(:,:,i-1);
%     A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)+1)*(M-1)+1:(2*(i)+1)*(M-1))=-eye(2*(M-1),2*(M-1));
% end
% 
% A(2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1),(2*(K-2)+1)*(M-1)+1:(2*(K-1)+1)*(M-1))=V(1:M-1,:,K-1);
% %A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
% A=sparse(A);
% YYY=lsqminnorm(-A,F',1.e-6);
% % YYY=-A\F';
% % Yy=Yy+YYY';
% 
% YYy=Yy+YYY';
% Yy=YYy;
% m=1;
% norm(F,inf);
% 
% while (norm(F,inf)/max(1,norm(Yy,inf))>10^(-4) && m<20)
%  
%   Y(:,1)=Yy(1:M-1);
% for i=2:2*K-1
%  Y(:,i)=Yy((M-1)*(i-1)+1:i*(M-1));
% end
% 
%     
%     
% YY(:,1)=zeros((M-1),1);
% AP=zeros((M-1),K-1);
% AV=zeros((M-1),K-1);
% [YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde1dmatrix(tt(1),tt(2),N,Y(:,1),pp0(:,1),M);
% 
% for i=1:K-1
%     CV(:,i)=Y(:,2*i+1);
%     CP(:,i)=Y(:,2*i);
% end 
% Tt=zeros(2,K-1);
% for i=1:K-1
%     Tt(:,i)=[tt(i+1),tt(i+2)];
% end 
% parfor i=1:K-1
%     [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde1dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
% end
%    
% for i=2:K
%     YY(:,2*i)=AV(:,i-1);
%     YY(:,2*i+1)=AP(:,i-1);
% end 
% 
% for i=1:K-1
%     %p 
% F((2*(i-1))*(M-1)+1:(2*(i-1)+1)*(M-1))= YY(:,2*i+1)-Y(:,2*i);
%     %v
% F((2*i-1)*(M-1)+1:(2*i)*(M-1))= YY(:,2*i)-Y(:,2*i+1);
% end
% F((2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1)))=YY(:,2*K+1)-p0(1:M-1,K+1);
% 
% 
% A(1:M-1,1:M-1)=VV(1:M-1,M:2*(M-1),1); % partial p1 \partial v0
% A(1:M-1,M:2*(M-1))=-eye(M-1,M-1); % partial p1 \partial v1
% A(M:2*(M-1),1:M-1)=VV(M:2*(M-1),M:2*(M-1),1); %\partial v1 \partial v0
% A(M:2*(M-1),2*(M-1)+1:3*(M-1))=-eye(M-1,M-1); %\partial v1 \partial p1
% for i=2:K-1
%     A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)-1)*(M-1)+1:(2*(i-1)+1)*(M-1))=V(:,:,i-1);
%     A(2*(i-1)*(M-1)+1:2*i*(M-1),(2*(i-1)+1)*(M-1)+1:(2*(i)+1)*(M-1))=-eye(2*(M-1),2*(M-1));
% end
% 
% A(2*(K-1)*(M-1)+1:(2*(K-1)+1)*(M-1),(2*(K-2)+1)*(M-1)+1:(2*(K-1)+1)*(M-1))=V(1:M-1,:,K-1);
% %A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
% A=sparse(A);
% YYY=lsqminnorm(-A,F',1.e-6);
% % YYY=-A\F';
% % Yy=Yy+YYY';
% 
% YYy=Yy+YYY';
% Yy=YYy;
% 
% norm(F, inf)
% m=m+1;
% toc(a)
% end
%        
%     
% end
% 
% 
% 
%  Y(:,1)=Yy(1:M-1);
% for i=2:2*K-1
%  Y(:,i)=Yy((M-1)*(i-1)+1:i*(M-1));
% end
% 
% vv0=(Y(1:M-1))';
% for i=1:K-1
%   % p
%   pp(:,i)=Y((2*i-1)*(M-1)+1:2*i*(M-1));
%   % S
%   VVV(:,i)=Y(2*i*(M-1)+1:(2*i+1)*(M-1));
% end   
% % SS;
% % pp;
% 
% for i=1:K-1
%    pmed(1:M-1,i)=Y(:,2*i);
%     end 
%    pmed(M,i)=(1-sum(Y(:,2*i))*dx)/(dx);
%  ppp(:,2:K)=pmed;  
%  ppp(:,1)=p0(:,1);
%  ppp(:,K+1)=p0(:,K+1);
%  
% 
