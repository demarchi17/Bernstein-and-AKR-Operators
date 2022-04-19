clear all
close all

%---------------------------------------------------------
% Copyright by Ana Maria Acu and Stefano De Marchi
% April 2022
%
%  This script computes the Aldaz-Kouncev-Render 
%  and Bernstein approximating operators
%  of a function in [0,1] and in [0,1]^2
%
% Reference
% A. Acu, S. De Marchi and I. Rasa
% Aldaz-Kounchev-Render operators and their approximation properties 
% (submitted April 2022)
%---------------------------------------------------------

 N=input('give me N=');
 j=input('give me j (>=1)');
%j=2; 
%N=4;    %%%%%%% Put N=3 for examples 4.2 and 4.3

M=50;   %nr. evaluation points

disp('1 dimensional example')

f1d=@(x) -2/pi*x.*cos(pi*x/2)+4/pi^2*sin(pi*x/2); %function to approximate
x=linspace(0,1,M);   %evaluation points
y=x;


for n=j:N


% Aldaz points

t(1:j)=0;
for k=j:n,
    t(k+1)= (factorial(k)/factorial(k-j)*factorial(n-j)/factorial(n))^(1/j);
end

% Bernstein points
b=0:1/n:1;


% Bernstein and Aldaz approximations in 1d

for i=1:length(x),
    for k=0:n,
       bern(k+1)=nchoosek(n,k)*x(i)^k*(1-x(i))^(n-k);
    end
    bbA(i)=f1d(t)*bern';
    bbB(i)=f1d(b)*bern';
end


ff=f1d(x);
eA(n)=norm(ff-bbA,inf)/norm(ff);
eB(n)=norm(ff-bbB,inf)/norm(ff);
end

figure
plot(x,f1d(x),'r-',  x,bbB,'k--',x,bbA,'b-.');
legend('function','B-Approx','AKR-approx');

figure
semilogy(1:length(eA),eA,'b--', 1:length(eB), eB, 'k-.')
legend('relerr 1d AKR','relerr 1d Bernstein')
grid


clear n
%-------------------------
%  bivariate example
%-------------------------
disp('2 dimensional example')
ExampleNr=input('Example 4.2 type 2, Example 4.3 type 3, Example 4.4 type 4 : ')



for n=j:N
 if ExampleNr==2,
    m=n+1;
    f2d=@(x,y) exp(x.^2.*y.^2)-1;       %Example 4.2
end
 if ExampleNr==3
    m=n+1;
    f2d=@(x,y) tan(pi/4*(2.^(x.^j)-1).*y.^(3*j)); %Example 4.3
 end
 if ExampleNr==4
    m=n;
   f2d=@(x,y) (1/pi)^4*(-64*pi*(y+x).*cos(pi/4*(y+x))+(256-16*pi^2*x.*y).*sin(pi/4*(x+y))+...
   64*pi*x.*cos(pi*x/4)+64*pi*y.*cos(pi/4*y)-256*(sin(pi/4*x)+sin(pi/4*y))); %Example 4.4
 end
 if ExampleNr<1 | ExampleNr>4
   error('2d example number not studied!')
 end
    
    
% Bernstein points
bbx=0:1/n:1;
bby=0:1/m:1;


% Aldaz points along x
aax(1:j)=0;
for k=j:n
    aax(k+1)= (factorial(k)/factorial(k-j)*factorial(n-j)/factorial(n))^(1/j);
end

aay(1:j)=0;
for k=j:m
    aay(k+1)= (factorial(k)/factorial(k-j)*factorial(m-j)/factorial(m))^(1/j);
end

 [bx,by]=meshgrid(bbx,bby);   %grid of Bernstein points
 [ax,ay]=meshgrid(aax,aay);   %grid of Aldaz points
 
 
 f2dB=f2d(bx,by);
 f2dA=f2d(ax,ay);
 %%%%%%
 w=x;
 %t=x;
 
 for i=1:length(w),
  for k=0:n,
     Bx(i,k+1)=nchoosek(n,k)*w(i)^k*(1-w(i))^(n-k);
    end
 end
for i=1:length(w),
  for k=0:m,
     By(i,k+1)=nchoosek(m,k)*w(i)^k*(1-w(i))^(m-k);
    end
end
 
%By=Bx;
 
 %Bx=bernsteinMatrix(n,w);

 [X,Y]=meshgrid(x);
 F=f2d(X,Y);
 %mesh(X,Y,F)
 
 BBernstein=By*f2dB*Bx';
 BAldaz=By*f2dA*Bx';
 
 errBer2d(n)=norm(BBernstein-F)/norm(F);
 errAKR2d(n)=norm(BAldaz-F)/norm(F);

end

% figure
% semilogy(1:length(errBer2d),errBer2d,'b-.',1:length(errAKR2d),errAKR2d,'r')
% legend('relerr Bernstein 2d','relerr AKR 2d')
% grid

subindex=[num2str(n),num2str(m),num2str(j)];
AKRstr=['B_{',subindex,'}f-f'];
Bstr=['B_{',num2str(n),num2str(m),'}f-f'];
BAKRstr=['B_{',subindex,'}f-','B_{',num2str(n),num2str(m),'}f'];

figure
surfc(X,Y,BAldaz-F)

title(AKRstr)
view(-60,30)

figure
surfc(X,Y,BBernstein-F)
title(Bstr)
view(-60,30)

figure
surfc(X,Y,BBernstein-BAldaz)
title(BAKRstr)
view(-60,30)
