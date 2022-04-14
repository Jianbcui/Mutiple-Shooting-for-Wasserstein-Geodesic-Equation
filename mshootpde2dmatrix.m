function [V,p,C]=mshootpde2dmatrix(ta,tb,N,V00,pp00,M)
% Neumann boundary condition for p
% Dirchlet boundary condition for v
dt=(tb-ta)/N;
t=[ta:dt:tb];

dx=0.2;

c=dt/dx; 
% regularization of density
alpha=0*10^-2;

% form density matrix 
P0(1:M^2-1)=pp00;
P0(M^2)=(1-(sum(pp00))*(dx)^2)/(dx)^2;
p=zeros(M+2,M+2);
p(2:M+1,2:M+1)=reshape(P0,[M,M]);
p(2:M+1,2:M+1)=p(2:M+1,2:M+1)';

%produce vx, vy, zero Dirchlet boundary condition for v
vx=zeros(M,M+1);
vy=zeros(M+1,M);
V=(reshape(V00,[M-1,M+1]))';
vx(:,2:M)=V(1:M,:);
vy(2:M,1)=V(M+1,:);
for j=2:M
vy(2:M,j)=vy(2:M,1)-sum(vx(1:M-1,1:j),2)+sum(vx(2:M,1:j),2);
end


% iteration number for fixed point iteration, w=1 is the Euler method  
w=1;
% change K get the density at different time 
% Average probability weight 
for l=1:N
    peuler=p;
    vxeuler=vx;
    vyeuler=vy;

    for m=1:w
        p(1,2:M+1)=p(2,2:M+1);
        p(M+2,2:M+1)=p(M+1,2:M+1);
        p(2:M+1,1)=p(2:M+1,2);
        p(2:M+1,M+2)=p(2:M+1,M+1);
        
 p(2:M+1,2:M+1)=peuler(2:M+1,2:M+1)+c*1/2*(p(2:M+1,1:M)+p(2:M+1,2:M+1)).*vx(:,1:M)-c*1/2*(p(2:M+1,2:M+1)+p(2:M+1,3:M+2)).*vx(:,2:M+1) ...
     +c*1/2*(p(1:M,2:M+1)+p(2:M+1,2:M+1)).*vy(1:M,:)-c*1/2*(p(2:M+1,2:M+1)+p(3:M+2,2:M+1)).*vy(2:M+1,:) ... 
     +alpha*(p(2:M+1,3:M+2)-2*p(2:M+1,2:M+1)+p(2:M+1,1:M)) ...
     +alpha*(p(1:M,2:M+1)-2*p(2:M+1,2:M+1)+p(3:M+2,2:M+1));    
 
               
 vy(2:M,:)=vyeuler(2:M,:)+1/4*c*dx*(vy(1:M-1,:)).^2-1/4*c*dx*(vy(3:M+1,:)).^2 ...
           +1/4*c*dx*(vx(1:M-1,1:M)).^2+1/4*c*dx*(vx(1:M-1,2:M+1)).^2 ...
           -1/4*c*dx*(vx(2:M,1:M)).^2-1/4*c*dx*(vx(2:M,2:M+1)).^2;
 
 vx(:,2:M)=vxeuler(:,2:M)+1/4*c*dx*(vx(1:M,1:M-1)).^2-1/4*c*dx*(vx(1:M,3:M+1)).^2 ...
           +1/4*c*dx*(vy(1:M,1:M-1)).^2-1/4*c*dx*(vy(1:M,2:M)).^2 ... 
           +1/4*c*dx*(vy(2:M+1,1:M-1)).^2-1/4*c*dx*(vy(2:M+1,2:M)).^2;
    end
      % p(find(p<=0))=10^-5;
end


%output for p
pend0=p(2:M+1,2:M+1)';
pend1=pend0(:);
pend=pend1(1:M^2-1);

%output for vx and vy
vend0=vx(:,2:M)';
vend=vend0(:);
vend(M^2-M+1:M^2-1)=vy(2:M,1);

% jacobi matrix 

% produce p^\epsilon, v
epsilon1=10^-6;
 XXx=repmat(pp00,1,M^2-1);
 YYy=repmat(V00,1,M^2-1);
 for l=1:(M^2-1)
 XXx(l,l)=XXx(l,l)+epsilon1;
 end 
 XXx(M^2,:)=(1-(sum(XXx(1:M^2-1,:),1))*(dx)^2)/(dx)^2;  
p=zeros(M+2,M+2,M^2-1);
% p^\epsilon 
p(2:M+1,2:M+1,:)=permute(reshape(XXx,[M,M,M^2-1]),[2,1,3]);
% v
vx=zeros(M,M+1,M^2-1);
vy=zeros(M+1,M,M^2-1);
V=permute(reshape(YYy,[M-1,M+1,M^2-1]),[2,1,3]);
vx(:,2:M,:)=V(1:M,:,:);
vy(2:M,1,:)=V(M+1,:,:);
for j=2:M
vy(2:M,j,:)=vy(2:M,1,:)-sum(vx(1:M-1,1:j,:),2)+sum(vx(2:M,1:j,:),2);
end


% iteration number for fixed point iteration, w=1 is the Euler method  
%w=1;
for i=1:N
     peuler=p;
    vxeuler=vx;
    vyeuler=vy;
    
    for m=1:w
        p(1,2:M+1,:)=p(2,2:M+1,:);
        p(M+2,2:M+1,:)=p(M+1,2:M+1,:);
        p(2:M+1,1,:)=p(2:M+1,2,:);
        p(2:M+1,M+2,:)=p(2:M+1,M+1,:);
 p(2:M+1,2:M+1,:)=peuler(2:M+1,2:M+1,:)+c*1/2*(p(2:M+1,1:M,:)+p(2:M+1,2:M+1,:)).*vx(:,1:M,:)-c*1/2*(p(2:M+1,2:M+1,:)+p(2:M+1,3:M+2,:)).*vx(:,2:M+1,:) ...
     +c*1/2*(p(1:M,2:M+1,:)+p(2:M+1,2:M+1,:)).*vy(1:M,:,:)-c*1/2*(p(2:M+1,2:M+1,:)+p(3:M+2,2:M+1,:)).*vy(2:M+1,:,:) ...
       +alpha*(p(2:M+1,3:M+2,:)-2*p(2:M+1,2:M+1,:)+p(2:M+1,1:M,:)) ...
     +alpha*(p(1:M,2:M+1,:)-2*p(2:M+1,2:M+1,:)+p(3:M+2,2:M+1,:));  

 vy(2:M,:,:)=vyeuler(2:M,:,:)+1/4*c*dx*(vy(1:M-1,:,:)).^2-1/4*c*dx*(vy(3:M+1,:,:)).^2 ...
           +1/4*c*dx*(vx(1:M-1,1:M,:)).^2+1/4*c*dx*(vx(1:M-1,2:M+1,:)).^2 ...
           -1/4*c*dx*(vx(2:M,1:M,:)).^2-1/4*c*dx*(vx(2:M,2:M+1,:)).^2;
        vx(:,2:M,:)=vxeuler(:,2:M,:)+1/4*c*dx*(vx(1:M,1:M-1,:)).^2-1/4*c*dx*(vx(1:M,3:M+1,:)).^2 ...
           +1/4*c*dx*(vy(1:M,1:M-1,:)).^2-1/4*c*dx*(vy(1:M,2:M,:)).^2 ... 
           +1/4*c*dx*(vy(2:M+1,1:M-1,:)).^2-1/4*c*dx*(vy(2:M+1,2:M,:)).^2;
   % p(find(p<0))=10^-5;
    end
end
 
% partial p \partial p_0, \partial v \partial p_0

ppend0=permute(p(2:M+1,2:M+1,:),[2,1,3]);
ppend1=reshape(ppend0,[],M^2-1);
ppend=ppend1(1:M^2-1,:);
vvend0=permute(vx(:,2:M,:),[2,1,3]);
vvend(1:M^2-M,:)=reshape(vvend0(:),[],M^2-1);
vvend(M^2-M+1:M^2-1,:)=vy(2:M,1,:);

Ppend=repmat(pend,1,M^2-1); 
Vvend=repmat(vend,1,M^2-1);  
  C1=(ppend-Ppend)/epsilon1;
  C2=(vvend-Vvend)/epsilon1;    

   for l=1:M^2-1
   C(1:M^2-1,l)=C1(:,l);
   C(M^2:2*M^2-2,l)=C2(:,l);
   end 
   
% jacobi matrix part 2

% produce p, v^epsilon
epsilon1=10^-6;
 XXx=repmat(pp00,1,M^2-1);
 YYy=repmat(V00,1,M^2-1);
 for l=1:(M^2-1)
 YYy(l,l)=YYy(l,l)+epsilon1;
 end 
 XXx(M^2,:)=(1-(sum(XXx(1:M^2-1,:),1))*(dx)^2)/(dx)^2; 
 % p
p=zeros(M+2,M+2,M^2-1);
p(2:M+1,2:M+1,:)=permute(reshape(XXx,[M,M,M^2-1]),[2,1,3]);
% v^epsilon
vx=zeros(M,M+1,M^2-1);
vy=zeros(M+1,M,M^2-1);
V=permute(reshape(YYy,[M-1,M+1,M^2-1]),[2,1,3]);
vx(:,2:M,:)=V(1:M,:,:);
vy(2:M,1,:)=V(M+1,:,:);
for j=2:M
vy(2:M,j,:)=vy(2:M,1,:)-sum(vx(1:M-1,1:j,:),2)+sum(vx(2:M,1:j,:),2);
end


% iteration number for fixed point iteration, w=1 is the Euler method  
%w=1;
for i=1:N
     peuler=p;
    vxeuler=vx;
    vyeuler=vy;
    for m=1:w
          p(1,2:M+1,:)=p(2,2:M+1,:);
        p(M+2,2:M+1,:)=p(M+1,2:M+1,:);
        p(2:M+1,1,:)=p(2:M+1,2,:);
        p(2:M+1,M+2,:)=p(2:M+1,M+1,:);    
 p(2:M+1,2:M+1,:)=peuler(2:M+1,2:M+1,:)+c*1/2*(p(2:M+1,1:M,:)+p(2:M+1,2:M+1,:)).*vx(:,1:M,:)-c*1/2*(p(2:M+1,2:M+1,:)+p(2:M+1,3:M+2,:)).*vx(:,2:M+1,:) ...
     +c*1/2*(p(1:M,2:M+1,:)+p(2:M+1,2:M+1,:)).*vy(1:M,:,:)-c*1/2*(p(2:M+1,2:M+1,:)+p(3:M+2,2:M+1,:)).*vy(2:M+1,:,:)...
     +alpha*(p(2:M+1,3:M+2,:)-2*p(2:M+1,2:M+1,:)+p(2:M+1,1:M,:)) ...
     +alpha*(p(1:M,2:M+1,:)-2*p(2:M+1,2:M+1,:)+p(3:M+2,2:M+1,:));  

 vy(2:M,:,:)=vyeuler(2:M,:,:)+1/4*c*dx*(vy(1:M-1,:,:)).^2-1/4*c*dx*(vy(3:M+1,:,:)).^2 ...
           +1/4*c*dx*(vx(1:M-1,1:M,:)).^2+1/4*c*dx*(vx(1:M-1,2:M+1,:)).^2 ...
           -1/4*c*dx*(vx(2:M,1:M,:)).^2-1/4*c*dx*(vx(2:M,2:M+1,:)).^2;
        vx(:,2:M,:)=vxeuler(:,2:M,:)+1/4*c*dx*(vx(1:M,1:M-1,:)).^2-1/4*c*dx*(vx(1:M,3:M+1,:)).^2 ...
           +1/4*c*dx*(vy(1:M,1:M-1,:)).^2-1/4*c*dx*(vy(1:M,2:M,:)).^2 ... 
           +1/4*c*dx*(vy(2:M+1,1:M-1,:)).^2-1/4*c*dx*(vy(2:M+1,2:M,:)).^2;           
   % p(find(p<0))=10^-5;
    end
end
 
% partial p \partial p_0, \partial v \partial p_0

ppend0=permute(p(2:M+1,2:M+1,:),[2,1,3]);
ppend1=reshape(ppend0,[],M^2-1);
ppend=ppend1(1:M^2-1,:);
vvend0=permute(vx(:,2:M,:),[2,1,3]);
vvend(1:M^2-M,:)=reshape(vvend0(:),[],M^2-1);
vvend(M^2-M+1:M^2-1,:)=vy(2:M,1,:);

Ppend=repmat(pend,1,M^2-1); 
Vvend=repmat(vend,1,M^2-1);  
  C3=(ppend-Ppend)/epsilon1;
  C4=(vvend-Vvend)/epsilon1;    
   
   for l=1:M^2-1
   C(1:M^2-1,M^2-1+l)=C3(:,l);
   C(M^2:2*M^2-2,M^2-1+l)=C4(:,l);
   end 
   
     
   V=vend(:,1); %output vx (M*(M-1)), vy ((M-1)*1)
   p=pend(:,1); % p(1:M^2-1)
     C; % numerical Jacobi matrix 
   

