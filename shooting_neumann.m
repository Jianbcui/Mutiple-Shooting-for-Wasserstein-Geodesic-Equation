clear all
 
a=tic;
T=1;
% small interval
K=10;
N=30;
tt=[0:T/K:T];
dx=0.2;
[x,y]=meshgrid(-1:dx:3,-1:dx:3);
M=length(x);
% newton iteration error
epsilon=10^(-10);
% shift of the density
r=0.001;


% % fake terminal distribution with lambda=0.1
% p0(:,:,K+1)=exp(-10*(y-0.6).^2-10*(x-0.4).^2)+r;
%      B=sum(sum(p0(:,:,K+1))*(dx)^2);
%      p0(:,:,K+1)=p0(:,:,K+1)/B;
%  % real initial distribution    
%     p0(:,:,1)=exp(-10*(y-0.5).^2-10*(x-0.3).^2)+r;
%      B=sum(sum(p0(:,:,1))*(dx)^2);
%      p0(:,:,1)=p0(:,:,1)/B;
% 
% % exact Gaussian interpolation as initial guess, it is the exact solution with lambda=0.1
%  for j=1:K-1
%       p0(:,:,j+1)=exp(-10*(y-0.5+j/K*(y-0.5-y+0.5-0.1)).^2-10*(x-0.3+j/K*(x-0.3-x+0.3-0.1)).^2)+r;
%       B=sum(sum(p0(:,:,j+1))*(dx)^2);
%       p0(:,:,j+1)=p0(:,:,j+1)/B;
%       
%  end

% different gaussian target 
% fake terminal distribution with lambda=0.1
% p0(:,:,K+1)=exp(-1/2*10*(y-1.5).^2-1/2*20*(x-1.3).^2)+r;
%      B=sum(sum(p0(:,:,K+1))*(dx)^2);
%      p0(:,:,K+1)=p0(:,:,K+1)/B;
% % real initial distribution      
% p0(:,:,1)=exp(-1/2*5*(y-0.5).^2-1/2*10*(x-0.3).^2)+r;
%      B=sum(sum(p0(:,:,1))*(dx)^2);
%      p0(:,:,1)=p0(:,:,1)/B;
% 
% % exact Gaussian interpolation as initial guess, it is the exact solution with lambda=0.1
% lambda=0.1;
%  for j=1:K-1
%       p0(:,:,j+1)=exp(-1/2*10*(y*(1/sqrt(2)-j/K*(1/sqrt(2)-1)*lambda)+j/K*(1/sqrt(2)*0.5-1.5)*lambda-1/sqrt(2)*0.5).^2 ...
%                   -1/2*20*(x*(1/sqrt(2)-j/K*(1/sqrt(2)-1)*lambda)+j/K*(1/sqrt(2)*0.3-1.3)*lambda-1/sqrt(2)*0.3).^2)+r;
%       B=sum(sum(p0(:,:,j+1))*(dx)^2);
%       p0(:,:,j+1)=p0(:,:,j+1)/B;
%       
%  end
% lambda=0.4;
% % fake terminal distribution with lambda=0.1
%  p0(:,:,K+1)=exp(-1/2*10*(y*(1/sqrt(2)-(1/sqrt(2)-1)*lambda)+(1/sqrt(2)*0.5-1.5)*lambda-1/sqrt(2)*0.5).^2 ...
%                   -1/2*20*(x*(1/sqrt(2)-(1/sqrt(2)-1)*lambda)+(1/sqrt(2)*0.3-1.3)*lambda-1/sqrt(2)*0.3).^2)+r;
%       B=sum(sum(p0(:,:,K+1))*(dx)^2);
%       p0(:,:,K+1)=p0(:,:,K+1)/B;
% % real initial distribution      
% p0(:,:,1)=exp(-1/2*5*(y-0.5).^2-1/2*10*(x-0.3).^2)+r;
%      B=sum(sum(p0(:,:,1))*(dx)^2);
%      p0(:,:,1)=p0(:,:,1)/B;

% lapacian distribution

% p0(:,:,1)= exp(-5*abs(x-0.5)-5*abs(y-0.6))+r;
%      B=sum(sum(p0(:,:,1))*(dx)^2);
%      p0(:,:,1)=p0(:,:,1)/B;

% 
%  p0(:,:,1)= ones(M,M);
%       B=sum(sum(p0(:,:,1))*(dx)^2);
%       p0(:,:,1)=p0(:,:,1)/B;
% % 
% p0(:,:,1)=exp(-50*(y-0.5).^2-50*(x-0.3).^2)+r;
%      B=sum(sum(p0(:,:,1))*(dx)^2);
%      p0(:,:,1)=p0(:,:,1)/B;
     
p0(:,:,1)=exp(-5*(y-0.5).^2-5*(x-0.3).^2)+r;
    B=sum(sum(p0(:,:,1))*(dx)^2);
    p0(:,:,1)=p0(:,:,1)/B;
% 
% % exact Gaussian interpolation as initial guess, it is the exact solution with lambda=0.1
%  for j=1:K-1
%       p0(:,:,j+1)=exp(-10*(y-0.5+j/K*(y-0.5-y+0.5-0.1)).^2-10*(x-0.3+j/K*(x-0.3-x+0.3-0.1)).^2)+r;
%       B=sum(sum(p0(:,:,j+1))*(dx)^2);
%       p0(:,:,j+1)=p0(:,:,j+1)/B;
%       
%  end
lambda=0.5;

% p0(:,:,K+1)=p0(:,:,1)+(p0(:,:,K+1)-p0(:,:,1))*lambda;

% p0(:,:,K+1)=exp(-50*(y-0.5-1*lambda).^2-50*(x-0.3-1*lambda).^2)+r;
%     B=sum(sum(p0(:,:,K+1))*(dx)^2);
%     p0(:,:,K+1)=p0(:,:,K+1)/B;

p0(:,:,K+1)=exp(-5*(y-0.5-lambda).^2-5*(x-0.3-lambda).^2)+r;
    B=sum(sum(p0(:,:,K+1))*(dx)^2);
    p0(:,:,K+1)=p0(:,:,K+1)/B;
% 
%  p0(:,:,1)= ((x+1).^2.*(3-x).^2+(y+1).^2.*(3-y).^2);
%        B=sum(sum(p0(:,:,1))*(dx)^2);
%        p0(:,:,1)=p0(:,:,1)/B;
% 
% 
% 
% % p0(:,:,1)=exp(-50*(y-0.5).^2-50*(x-0.3).^2)+r;
% %      B=sum(sum(p0(:,:,1))*(dx)^2);
% %      p0(:,:,1)=p0(:,:,1)/B;
% % 
% % % exact Gaussian interpolation as initial guess, it is the exact solution with lambda=0.1
% %  for j=1:K-1
% %       p0(:,:,j+1)=exp(-10*(y-0.5+j/K*(y-0.5-y+0.5-0.1)).^2-10*(x-0.3+j/K*(x-0.3-x+0.3-0.1)).^2)+r;
% %       B=sum(sum(p0(:,:,j+1))*(dx)^2);
% %       p0(:,:,j+1)=p0(:,:,j+1)/B;
% %       
% %  end
% lambda=0.1;
% 
% p0(:,:,K+1)= exp(-10*abs(x-1)-10*abs(y-1))+r;
%      B=sum(sum(p0(:,:,K+1))*(dx)^2);
%      p0(:,:,K+1)=p0(:,:,K+1)/B;
     
% 


% 
% p0(:,:,K+1)=ones(M,M);
%     B=sum(sum(p0(:,:,K+1))*(dx)^2);
%     p0(:,:,K+1)=p0(:,:,K+1)/B;

% p0(:,:,K+1)=exp(-50*(y-0.5-1*lambda).^2-50*(x-0.3-1*lambda).^2)+r;
%     B=sum(sum(p0(:,:,K+1))*(dx)^2);
%     p0(:,:,K+1)=p0(:,:,K+1)/B;

 
% for j=1:K-1
%     p0(:,:,j+1)=p0(:,:,1)+j*(p0(:,:,K+1)-p0(:,:,1))/K*lambda;
% end
 for j=1:K-1
p0(:,:,j+1)=exp(-5*(y-0.5-lambda*j/K).^2-5*(x-0.5-lambda*j/K).^2)+r;
    B=sum(sum(p0(:,:,j+1))*(dx)^2);
    p0(:,:,j+1)=p0(:,:,j+1)/B;
    
 end
    
% p0(:,:,K+1)=p0(:,:,1)+K*(p0(:,:,K+1)-p0(:,:,1))/K*lambda;
%      
%      for j=1:K-1
%     p0(:,:,j+1)=p0(:,:,1)+j*(p0(:,:,K+1)-p0(:,:,1))/K;
%     end
%     
% p0(:,:,K+1)=p0(:,:,1)+K*(p0(:,:,K+1)-p0(:,:,1))/K;

% the exact constan velocity as initial guess, it is the exact solution of
% the fake terminal distribution with lambda=0.1
% for l=1:K 
% for j=1:M
% V0((M-1)*(j-1)+1:(M-1)*j,l)=1.5*ones(M-1,1)/K-1/sqrt(2)*0.5*ones(M-1,1)/K+0.1*x(1:M-1,1)*(1/sqrt(2)-1);
% end
% for k=1:M-1
% V0((M-1)*M+k,l)=(1.3-1/sqrt(2)*0.3)/K+0.1*y(k,1)*(1/sqrt(2)-1);
% end
% end 

for l=1:K 
for j=1:M
V0((M-1)*(j-1)+1:(M-1)*j,l)=0*ones(M-1,1)/K;
end
for k=1:M-1
V0((M-1)*M+k,l)=0*1/K;
end
end 


% initial guess for density and velocity 
X0(1:M)=p0(1,:,1);
for j=2:M
   X0((j-1)*M+1:j*M)=p0(j,:,1); 
end 
pp0(:,1)=X0;
for i=2:K+1
   X0(1:M)=p0(1,:,i);
   for j=2:M
   X0((j-1)*M+1:j*M)=p0(j,:,i); 
   end 
   pp0(:,i)=X0;
end 
%  V0=ones(M^2-1,K)/K;

% input vectors is $p,v$ at each subintervel
Y(:,1)=V0(:,1);
for i=1:K-1
  % p
  Y(:,2*i)=pp0(1:M^2-1,i+1);
  % v
  Y(:,2*i+1)=V0(:,i+1);
end 

% transform the input vectors into a big vector
 Yy(1:M^2-1)=Y(:,1);
for i=2:2*K-1
 Yy((M^2-1)*(i-1)+1:i*(M^2-1))=Y(:,i);
end


YY(:,1)=zeros((M+1)*(M-1),1);
% V1,P1 via V0, P0 in [0,dt]
AP=zeros((M+1)*(M-1),K-1);
AV=zeros((M+1)*(M-1),K-1);

 % VV=zeros((M+1)*(M-1),(M+1)*(M-1),K);
 %is the jacobi matrix of the euler
% method in first subinterval 
[YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde2dmatrix(tt(1),tt(2),N,V0(:,1),pp0(1:M^2-1,1),M);
% parfor is used for computing each subintervals at the same time 
for i=1:K-1
    CV(:,i)=Y(:,2*i+1);
    CP(:,i)=Y(:,2*i);
end 


Tt=zeros(2,K-1);
for i=1:K-1
    Tt(:,i)=[tt(i+1),tt(i+2)];
end 

% V=zeros((M+1)*(M-1),(M+1)*(M-1),K) is the jacobi matrix of the euler
% method in each subinterval 
parfor i=1:K-1
    [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde2dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
end

% Put output in each subinterval into YY
for i=2:K
    YY(:,2*i)=AV(:,i-1);
    YY(:,2*i+1)=AP(:,i-1);
end 

% perform newton method
% goal [V0 p1 V1 ... pN-1,VN-1]
% objective function 
for i=1:K-1
    %p 
F((2*(i-1))*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))= YY(:,2*i+1)-Y(:,2*i);
    %v
F((2*i-1)*(M^2-1)+1:(2*i)*(M^2-1))= YY(:,2*i)-Y(:,2*i+1);
end
F((2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1)))=YY(:,2*K+1)-pp0(1:M^2-1,K+1);
% F_{\lambda}=F/


% Form the jacobi matrix of multiple shooting 
A(1:M^2-1,1:M^2-1)=VV(1:M^2-1,M^2:2*(M^2-1),1); % partial p1 \partial v0
A(1:M^2-1,M^2:2*(M^2-1))=-eye(M^2-1,M^2-1); % partial p1 \partial v1
A(M^2:2*(M^2-1),1:M^2-1)=VV(M^2:2*(M^2-1),M^2:2*(M^2-1),1); %\partial v1 \partial v0
A(M^2:2*(M^2-1),2*(M^2-1)+1:3*(M^2-1))=-eye(M^2-1,M^2-1); %\partial v1 \partial p1
for i=2:K-1
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)-1)*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))=V(:,:,i-1);
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)+1)*(M^2-1)+1:(2*(i)+1)*(M^2-1))=-eye(2*(M^2-1),2*(M^2-1));
end

A(2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1),(2*(K-2)+1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1))=V(1:M^2-1,:,K-1);
%A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
A=sparse(A);
YYY=lsqminnorm(-A,F',1.e-6);
% YYY=-A\F';
% Yy=Yy+YYY';

YYy=Yy+YYY';

for i=1:K-1
  minp(1:M^2-1,i+1)=YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
end   

kk=0;
% while (min(min(minp(1:M^2-1,2:K)))<=0 && kk<10)
%   
%     
% for i=1:K-1
%   minp(1:M^2-1,i+1)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))+1/(2^kk)*(YYY((2*i-1)*(M^2-1)+1:2*i*(M^2-1)))';
% end   
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end
%   kk=1+kk;
% end


% 
% % minp(find(minp<=0))=10^-5;
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end


Yy=YYy;

% errp=ones(1,K);
% 
% for i=1:K-1
%     %p 
% errp(i)=max(abs(F((2*(i-1))*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))));
%     %v
% end
% errp(K)=max(abs(F((2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1)))));
% 
% max(errp)
% 
%  modY(:,1)=YYy(1:M^2-1);
% for i=2:2*K-1
%  modY(:,i)=YYy((M^2-1)*(i-1)+1:i*(M^2-1));
% end
% 
% 
% vv0=(Yy(1:M^2-1))';
% for i=1:K-1
%   % p
%   pp(:,i)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
%   % S
%   VVV(:,i)=Yy(2*i*(M^2-1)+1:(2*i+1)*(M^2-1));
% end   
% 
% 
% for i=1:K-1
%     for j=1:M-1
%    pmed(j,:,i)=modY(M*(j-1)+1:M*j,2*i);
%     end 
%    pmed(M,1:M-1,i)=modY(M*(M-1)+1:M*(M-1)+M-1,2*i);
%    pmed(M,M,i)=(1-sum(modY(:,2*i))*dx^2)/(dx)^2;
% end 

% 
%  modY(:,1)=YYy(1:M^2-1);
% for i=2:2*K-1
%  modY(:,i)=YYy((M^2-1)*(i-1)+1:i*(M^2-1));
% end

% e=max(max(abs(modY-Y)))
m=1
norm(F,inf);
epsilon2=0.002*10^(-1);


while (norm(F, inf)>epsilon2 && norm(F,inf)/norm(Yy,inf)>10^(-5) && m<10)

  Y(:,1)=Yy(1:M^2-1);
for i=2:2*K-1
 Y(:,i)=Yy((M^2-1)*(i-1)+1:i*(M^2-1));
end

 
YY(:,1)=zeros((M+1)*(M-1),1);
% V1,P1 via V0, P0 in [0,dt]
AP=zeros((M+1)*(M-1),K-1);
AV=zeros((M+1)*(M-1),K-1);
% VV=zeros((M+1)*(M-1),(M+1)*(M-1),K);
[YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde2dmatrix(tt(1),tt(2),N,Y(:,1),pp0(1:M^2-1,1),M);

for i=1:K-1
    CV(:,i)=Y(:,2*i+1);
    CP(:,i)=Y(:,2*i);
end 


Tt=zeros(2,K-1);
for i=1:K-1
    Tt(:,i)=[tt(i+1),tt(i+2)];
end 

parfor i=1:K-1
    [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde2dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
end

% V=zeros((M+1)*(M-1),(M+1)*(M-1),K) is the jacobi matrix of the euler
% method in each subinterval 

% updating YY 
for i=2:K
    YY(:,2*i)=AV(:,i-1);
    YY(:,2*i+1)=AP(:,i-1);
end 


% perform newton method
% goal [v0 p1 v1 ... pN-1,vN-1]
% objective function 
for i=1:K-1
    %p 
F((2*(i-1))*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))= YY(:,2*i+1)-Y(:,2*i);
    %v
F((2*i-1)*(M^2-1)+1:(2*i)*(M^2-1))= YY(:,2*i)-Y(:,2*i+1);
end
F((2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1)))=YY(:,2*K+1)-pp0(1:M^2-1,K+1);

% F/delta lambda



% jacobi matrix


A(1:M^2-1,1:M^2-1)=VV(1:M^2-1,M^2:2*(M^2-1),1); % partial p1 \partial v0
A(1:M^2-1,M^2:2*(M^2-1))=-eye(M^2-1,M^2-1); % partial p1 \partial v1
A(M^2:2*(M^2-1),1:M^2-1)=VV(M^2:2*(M^2-1),M^2:2*(M^2-1),1); %\partial v1 \partial v0
A(M^2:2*(M^2-1),2*(M^2-1)+1:3*(M^2-1))=-eye(M^2-1,M^2-1); %\partial v1 \partial p1
for i=2:K-1
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)-1)*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))=V(:,:,i-1);
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)+1)*(M^2-1)+1:(2*(i)+1)*(M^2-1))=-eye(2*(M^2-1),2*(M^2-1));
end
A(2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1),(2*(K-2)+1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1))=V(1:M^2-1,:,K-1);

A=sparse(A);
% YYY=-A\F';
%de=det(full(A))
YYY=lsqminnorm(-A,F',1.e-6);
%Yy=Yy+YYY';

YYy=Yy+YYY';

for i=1:K-1
  minp(1:M^2-1,i+1)=YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
end   

kk=0;
% while (min(min(minp(1:M^2-1,2:K)))<=0 && kk<10)
%   
%     
% for i=1:K-1
%   minp(1:M^2-1,i+1)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))+1/(2^kk)*(YYY((2*i-1)*(M^2-1)+1:2*i*(M^2-1)))';
% end   
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end
%   kk=1+kk;
% end
% while (min(min(minp(1:M^2-1,2:K)))<=0 && kk<10)
%   
%     
% for i=1:K-1
%   minp(1:M^2-1,i+1)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))+1/(2^kk)*(YYY((2*i-1)*(M^2-1)+1:2*i*(M^2-1)))';
% end   
% 
% YYy=Yy+1/(2^kk)*YYY';
% 
%   kk=1+kk;
% end
% 
% % minp(find(minp<=0))=10^-5;
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end
% minp(find(minp<=0))=10^-12;
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end
% 
%  modY(:,1)=YYy(1:M^2-1);
% for i=2:2*K-1
%  modY(:,i)=YYy((M^2-1)*(i-1)+1:i*(M^2-1));
% end
% 
% % 
% e=max(max(abs(modY-Y)))
% 

% 
% 
% vv0=(Yy(1:M^2-1))';
% for i=1:K-1
%   % p
%   pp(:,i)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
%   % S
%   VVV(:,i)=Yy(2*i*(M^2-1)+1:(2*i+1)*(M^2-1));
% end   
% 
% 
% for i=1:K-1
%     for j=1:M-1
%    pmed(j,:,i)=modY(M*(j-1)+1:M*j,2*i);
%     end 
%    pmed(M,1:M-1,i)=modY(M*(M-1)+1:M*(M-1)+M-1,2*i);
%    pmed(M,M,i)=(1-sum(modY(:,2*i))*dx^2)/(dx)^2;
% end 






Yy=YYy;

norm(F, inf)
m=m+1;
toc(a)
end

% first step finished 
    lambda=0.5;
% second step repeat the initial guess 
for q=1:2
    lambda=lambda+0.5 /2
    
    % using the date in last step as the new initial guess    
V0(:,1)=(Yy(1:M^2-1))';
for i=1:K-1
  % p
  pp0(1:M^2-1,i+1)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
  % S
  V0(:,i+1)=Yy(2*i*(M^2-1)+1:(2*i+1)*(M^2-1));
end   

Y(:,1)=V0(:,1);
for i=1:K-1
  % p
  Y(:,2*i)=pp0(1:M^2-1,i+1);
  % v
  Y(:,2*i+1)=V0(:,i+1);
end 

% new fake target distribution 
%  p0(:,:,K+1)=exp(-1/2*10*(y*(1/sqrt(2)-(1/sqrt(2)-1)*lambda)+(1/sqrt(2)*0.5-1.5)*lambda-1/sqrt(2)*0.5).^2 ...
%                   -1/2*20*(x*(1/sqrt(2)-(1/sqrt(2)-1)*lambda)+(1/sqrt(2)*0.3-1.3)*lambda-1/sqrt(2)*0.3).^2)+r;
%       B=sum(sum(p0(:,:,K+1))*(dx)^2);
%       p0(:,:,K+1)=p0(:,:,K+1)/B;
% 
% p0(:,:,K+1)=exp(-50*(y-1.5+lambda).^2-50*(x-1.3+lambda).^2)+r;
%     B=sum(sum(p0(:,:,K+1))*(dx)^2);
%     p0(:,:,K+1)=p0(:,:,K+1)/B;
% 
p0(:,:,j+1)=exp(-5*(y-0.5-lambda).^2-5*(x-0.3-lambda).^2)+r;
    B=sum(sum(p0(:,:,K+1))*(dx)^2);
    p0(:,:,K+1)=p0(:,:,K+1)/B;
%  p0(:,:,K+1)=p0(:,:,1)+(p0(:,:,K+1)-p0(:,:,1))*lambda;


% p0(:,:,K+1)=p0(:,:,1)+(p0(:,:,K+1)-p0(:,:,1))*lambda;
%      B=sum(sum(p0(:,:,K+1))*(dx)^2);
%      p0(:,:,K+1)=p0(:,:,K+1)/B;
% p0(:,:,K+1)= exp(-10*abs(x-1*lambda)-10*abs(y-1*lambda))+r;
%      B=sum(sum(p0(:,:,K+1))*(dx)^2);
%      p0(:,:,K+1)=p0(:,:,K+1)/B;

% p0(:,:,K+1)=exp(-50*(y-1.5+lambda).^2-50*(x-1.3+lambda).^2)+r;
%     B=sum(sum(p0(:,:,K+1))*(dx)^2);
%     p0(:,:,K+1)=p0(:,:,K+1)/B;
%     p0(:,:,K+1)=exp(-50*(y-1.5+lambda).^2-50*(x-1.3+lambda).^2)+r;
%     B=sum(sum(p0(:,:,K+1))*(dx)^2);
%     p0(:,:,K+1)=p0(:,:,K+1)/B;
    
%      
% 
 %p0(:,:,K+1)=p0(:,:,1)+(p0(:,:,K+1)-p0(:,:,1))*lambda;

    
    
X0(1:M)=p0(1,:,K+1);
for j=2:M
   X0((j-1)*M+1:j*M)=p0(j,:,K+1); 
end 
pp0(:,K+1)=X0;

% new initial guess 


% the old value of [V0,p1,V1,\cdots,pN-1,VN-1]

 Yy(1:M^2-1)=Y(:,1);
for i=2:2*K-1
 Yy((M^2-1)*(i-1)+1:i*(M^2-1))=Y(:,i);
end


% form the new value of [V0,p1,V1,\cdots,pN-1,VN-1] by forward euler 
YY(:,1)=zeros((M+1)*(M-1),1);
% S1,P1 via S0, P0 in [0,dt]
AP=zeros((M+1)*(M-1),K-1);
AV=zeros((M+1)*(M-1),K-1);
% V=zeros((M+1)*(M-1),(M+1)*(M-1),K);
[YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde2dmatrix(tt(1),tt(2),N,V0(:,1),pp0(1:M^2-1,1),M);

% parfor computations 

for i=1:K-1
    CV(:,i)=Y(:,2*i+1);
    CP(:,i)=Y(:,2*i);
end 
% parfor time

Tt=zeros(2,K-1);
for i=1:K-1
    Tt(:,i)=[tt(i+1),tt(i+2)];
end 

parfor i=1:K-1
    [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde2dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
end

for i=2:K
    YY(:,2*i)=AV(:,i-1);
    YY(:,2*i+1)=AP(:,i-1);
end 

% perform newton method
% goal [V0 p1 V1 ... pN-1,VN-1]

% objective function 
for i=1:K-1
    %p 
F((2*(i-1))*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))= YY(:,2*i+1)-Y(:,2*i);
    %v
F((2*i-1)*(M^2-1)+1:(2*i)*(M^2-1))= YY(:,2*i)-Y(:,2*i+1);
end
F((2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1)))=YY(:,2*K+1)-pp0(1:M^2-1,K+1);

% jacobi matrix

A(1:M^2-1,1:M^2-1)=VV(1:M^2-1,M^2:2*(M^2-1),1); % partial p1 \partial v0
A(1:M^2-1,M^2:2*(M^2-1))=-eye(M^2-1,M^2-1); % partial p1 \partial v1
A(M^2:2*(M^2-1),1:M^2-1)=VV(M^2:2*(M^2-1),M^2:2*(M^2-1),1); %\partial v1 \partial v0
A(M^2:2*(M^2-1),2*(M^2-1)+1:3*(M^2-1))=-eye(M^2-1,M^2-1); %\partial v1 \partial p1
for i=2:K-1
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)-1)*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))=V(:,:,i-1);
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)+1)*(M^2-1)+1:(2*(i)+1)*(M^2-1))=-eye(2*(M^2-1),2*(M^2-1));
end

A(2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1),(2*(K-2)+1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1))=V(1:M^2-1,:,K-1);
%A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
A=sparse(A);
% YYY=-A\F';
% Yy=Yy+YYY';
%de=det(full(A))
YYY=lsqminnorm(-A,F',1.e-6);
YYy=Yy+YYY';
kk=0;

Yy=YYy;
m=1;
norm(F,inf)
epsilon2=0.002*10^(-1);

% continue iteration until convergence 
while (norm(F, inf)>epsilon2 && norm(F,inf)/norm(Yy,inf)>10^(-3) && m<10)% stopping critera

  Y(:,1)=Yy(1:M^2-1);
for i=2:2*K-1
 Y(:,i)=Yy((M^2-1)*(i-1)+1:i*(M^2-1));
end


 
YY(:,1)=zeros((M+1)*(M-1),1);
% S1,P1 via S0, P0 in [0,dt]
AP=zeros((M+1)*(M-1),K-1);
AV=zeros((M+1)*(M-1),K-1);
% V=zeros((M+1)*(M-1),(M+1)*(M-1),K);
[YY(:,2),YY(:,3),VV(:,:,1)]=mshootpde2dmatrix(tt(1),tt(2),N,Y(:,1),pp0(1:M^2-1,1),M);
% parfor time
for i=1:K-1
    CV(:,i)=Y(:,2*i+1);
    CP(:,i)=Y(:,2*i);
end 


Tt=zeros(2,K-1);
for i=1:K-1
    Tt(:,i)=[tt(i+1),tt(i+2)];
end 

parfor i=1:K-1
    [AV(:,i),AP(:,i),V(:,:,i)]=mshootpde2dmatrix(Tt(1,i),Tt(2,i),N,CV(:,i),CP(:,i),M);
end

for i=2:K
    YY(:,2*i)=AV(:,i-1);
    YY(:,2*i+1)=AP(:,i-1);
end 

% perform newton method
% goal [V0 p1 V1 ... pN-1,VN-1]
% objective function 
for i=1:K-1
    %p 
F((2*(i-1))*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))= YY(:,2*i+1)-Y(:,2*i);
    %v
F((2*i-1)*(M^2-1)+1:(2*i)*(M^2-1))= YY(:,2*i)-Y(:,2*i+1);
end
F((2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1)))=YY(:,2*K+1)-pp0(1:M^2-1,K+1);
% jacobi matrix
A(1:M^2-1,1:M^2-1)=VV(1:M^2-1,M^2:2*(M^2-1),1); % partial p1 \partial v0
A(1:M^2-1,M^2:2*(M^2-1))=-eye(M^2-1,M^2-1); % partial p1 \partial v1
A(M^2:2*(M^2-1),1:M^2-1)=VV(M^2:2*(M^2-1),M^2:2*(M^2-1),1); %\partial v1 \partial v0
A(M^2:2*(M^2-1),2*(M^2-1)+1:3*(M^2-1))=-eye(M^2-1,M^2-1); %\partial v1 \partial p1
for i=2:K-1
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)-1)*(M^2-1)+1:(2*(i-1)+1)*(M^2-1))=V(:,:,i-1);
    A(2*(i-1)*(M^2-1)+1:2*i*(M^2-1),(2*(i-1)+1)*(M^2-1)+1:(2*(i)+1)*(M^2-1))=-eye(2*(M^2-1),2*(M^2-1));
end
A(2*(K-1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1),(2*(K-2)+1)*(M^2-1)+1:(2*(K-1)+1)*(M^2-1))=V(1:M^2-1,:,K-1);
%A(2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1),2*(K-1)*(M+1)+1:(2*(K-1)+1)*(M+1))=-eye(M+1,M+1)
A=sparse(A);
% YYY=-A\F';
% Yy=Yy+YYY';

%de=det(full(A))
YYY=lsqminnorm(-A,F',1.e-6);
 YYy=Yy+YYY';


% while (min(min(minp(1:M^2-1,2:K)))<=0 && kk<10)
%   
%     
% for i=1:K-1
%   minp(1:M^2-1,i+1)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))+1/(2^kk)*(YYY((2*i-1)*(M^2-1)+1:2*i*(M^2-1)))';
% end   
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end
%   kk=1+kk;
% end
for i=1:K-1
  minp(1:M^2-1,i+1)=YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
end   


% minp(find(minp<=0))=10^-5;



% minp(find(minp<=0))=10^-5;

for i=1:K-1
YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
end
% minp(find(minp<=0))=10^-12;
% 
% for i=1:K-1
% YYy((2*i-1)*(M^2-1)+1:2*i*(M^2-1))=minp(1:M^2-1,i+1);
% end

%  modY(:,1)=YYy(1:M^2-1);
% for i=2:2*K-1
%  modY(:,i)=YYy((M^2-1)*(i-1)+1:i*(M^2-1));
% end

% 
% e=max(max(abs(modY-Y)))

Yy=YYy;

 Y(:,1)=Yy(1:M^2-1);
for i=2:2*K-1
 Y(:,i)=Yy((M^2-1)*(i-1)+1:i*(M^2-1));
end


vv0=(Yy(1:M^2-1))';
for i=1:K-1
  % p
  pp(:,i)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
  % S
  VVV(:,i)=Yy(2*i*(M^2-1)+1:(2*i+1)*(M^2-1));
end   


for i=1:K-1
    for j=1:M-1
   pmed(j,:,i)=Y(M*(j-1)+1:M*j,2*i);
    end 
   pmed(M,1:M-1,i)=Y(M*(M-1)+1:M*(M-1)+M-1,2*i);
   pmed(M,M,i)=(1-sum(Y(:,2*i))*dx^2)/(dx)^2;
end 



mp=min(min(min(pmed))) % lower bound of the density 
norm(F, inf)
m=m+1;
toc(a)
end

end

% final output

% for j=1:K-1
%      p0(:,:,K+1)=exp(-1/2*10*(y*(1/sqrt(2)-(1/sqrt(2)-1)*lambda)+(1/sqrt(2)*0.5-1.5)*lambda-1/sqrt(2)*0.5).^2 ...
%                   -1/2*20*(x*(1/sqrt(2)-(1/sqrt(2)-1)*lambda)+(1/sqrt(2)*0.3-1.3)*lambda-1/sqrt(2)*0.3).^2)+r;
%       B=sum(sum(p0(:,:,K+1))*(dx)^2);
%       p0(:,:,K+1)=p0(:,:,K+1)/B;
%       
% end

 Y(:,1)=Yy(1:M^2-1);
for i=2:2*K-1
 Y(:,i)=Yy((M^2-1)*(i-1)+1:i*(M^2-1));
end

vv0=(Yy(1:M^2-1))';
for i=1:K-1
  % p
  pp(:,i)=Yy((2*i-1)*(M^2-1)+1:2*i*(M^2-1));
  % S
  VVV(:,i)=Yy(2*i*(M^2-1)+1:(2*i+1)*(M^2-1));
end   



for i=1:K-1
    for j=1:M-1
   pmed(j,:,i)=Y(M*(j-1)+1:M*j,2*i);
    end 
   pmed(M,1:M-1,i)=Y(M*(M-1)+1:M*(M-1)+M-1,2*i);
   pmed(M,M,i)=(1-sum(Y(:,2*i))*dx^2)/(dx)^2;
end 
