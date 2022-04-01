clear all
close all

%---------------------------------------------------------
% Copyright by Ana Maria Acu and Stefano De Marchi
% March 2022
%
%  This script compute the Aldaz and Bernstein approximation
%  of a function on a triangle
%---------------------------------------------------------

% n=input('give me n=');
% j=input('give me j (>=1)');
j=2;
m=100; N=15;
x=linspace(0,1,m);   %evaluation points
y=x;

%f2d=@(x,y) exp(x.^2+y.^2);
f2d=@(x,y) exp(-x-y);

for n=j:N
% Bernstein points 1d
b=0:1/n:1;
% Aldaz points 1d
a(1:j)=0;
for k=j:n
    a(k+1)= (factorial(k)/factorial(k-j)*factorial(n-j)/factorial(n))^(1/j);
end

 [bx,by]=meshgrid(b);   %grid of Bernstein points
 [ax,ay]=meshgrid(a);   %grid of Aldaz points

% %-------------- Bernstein and AKR on the unit right triangle ----------%

btx=triu(fliplr(bx)); bty=triu(by); %B-points on triangle
atx=triu(fliplr(ax)); aty=triu(ay); %AKR points on triangle

%  %--- Function values at B and AKR points of the T
f2dBT=f2d(btx,bty);
f2dAT=f2d(atx,aty);

%  evaluation points on the T
[X,Y]=meshgrid(x);
ex=triu(fliplr(X)); ey=triu(Y);
ex=fliplr(ex); ey=fliplr(ey);

FT=f2d(ex,ey);  %function at evaluation on T


%x=fliplr(ex(1,:)); y=ey(:,m);

BBxyT=zeros(m,m);
BAKRxyT=zeros(m,m);
for r=1:m
for s=1:m-r+1
   for i=0:n
   ia=sqrt(i*(i-1)/(n*(n-1)));
   for l=0:n-i
     ff=factorial(n)/(factorial(i)*factorial(l)*factorial(n-i-l));    
     FB(i+1,l+1)=ff*f2d(i/n,l/n);
     CBil=FB(i+1,l+1);
     ja=sqrt(l*(l-1)/(n*(n-1)));
     FA(i+1,l+1)=ff*f2d(ia,ja);
     CAil=FA(i+1,l+1);      

%      BBxyT(r,s)=CBil*ex(r,s)^i*ey(r,s)^l*(1-ex(r,s)-ey(r,s))^(n-i-l);
%      BAKRxyT(r,s)=CAil*ex(r,s)^i*ey(r,s)^l*(1-ex(r,s)-ey(r,s))^(n-i-l);
     xy=x(r)^i*y(s)^l*(1-x(r)-y(s))^(n-i-l);
     BBxyT(r,s)=BBxyT(r,s)+CBil*xy;
     BAKRxyT(r,s)=BAKRxyT(r,s)+CAil*xy;
    end
  end
 end
end

FTT=triu(fliplr(FT));
BB=fliplr(BBxyT); BAKR=fliplr(BAKRxyT);
errBer2dT(n)=norm(BB-FTT)/norm(FTT);
errAKR2dT(n)=norm(BAKR-FTT)/norm(FTT);

end

subplot(3,1,1)
%surf(x,y,FTT)
s = surf(x,y,FTT,'FaceAlpha',0.7);
s.EdgeColor='none';
title('exp(-x-y)')

subplot(3,1,2)
%surf(x,y,BB)
s = surf(x,y,BB,'FaceAlpha',0.7);
s.EdgeColor='none';
title('B-approximant deg 15')

subplot(3,1,3)
%surf(x,y,BAKR)
s = surf(x,y,BAKR,'FaceAlpha',0.7);
s.EdgeColor='none';
title('AKR-approximant deg 15')


figure
plot(btx,bty,'o',atx,aty,'*')
legend('B-points','AKR points for j=2 on the triangle')

figure
semilogy(1:length(errBer2dT),errBer2dT,'b-.',1:length(errAKR2dT),errAKR2dT,'r')
legend('relerr Bernstein Triangle','relerr AKR Triangle')
grid
