clear all
close all

%---------------------------------------------------------
% Copyright by Ana Maria Acu and Stefano De Marchi
% March 2022
%
%  This script compute the Aldaz and Bernstein approximation
%  of a function in 1d and 2d as well as cubature
%---------------------------------------------------------

% n=input('give me n=');
% j=input('give me j (>=1)');
j=2;
m=100; N=20;
f1d=@(x) -2/pi*x.*cos(pi*x/2)+4/pi^2*sin(pi*x/2); %function to approximate
x=linspace(0,1,m);   %evaluation points
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


plot(x,f1d(x),'r', x,bbA,'b', x,bbB,'k');

ff=f1d(x);
eA(n)=norm(ff-bbA,inf)/norm(ff);
eB(n)=norm(ff-bbB,inf)/norm(ff);
end

semilogy(1:length(eA),eA,'r', 1:length(eB), eB, 'b-.')
legend('relerr 1d AKR','relerr 1d Bernstein')
grid


clear n
%
%  bivariate example
%
f2d=@(x,y) exp(x.^2.*y.^2);
%f2d=@(x,y) exp(-x-y);

for n=j:N
% figure(2)
% Bernstein points
b=0:1/n:1;


% Aldaz points
a(1:j)=0;
for k=j:n
    a(k+1)= (factorial(k)/factorial(k-j)*factorial(n-j)/factorial(n))^(1/j);
end

 [bx,by]=meshgrid(b);   %grid of Bernstein points
 [ax,ay]=meshgrid(a);   %grid of Aldaz points
 
 
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
By=Bx;
 
 %Bx=bernsteinMatrix(n,w);

 [X,Y]=meshgrid(x);
 F=f2d(X,Y);
 %mesh(X,Y,F)
 
 BBernstein=Bx*f2dB*By';
 BAldaz=Bx*f2dA*By';
 
 errBer2d(n)=norm(BBernstein-F)/norm(F);
 errAKR2d(n)=norm(BAldaz-F)/norm(F);



end

figure
semilogy(1:length(errBer2d),errBer2d,'b-.',1:length(errAKR2d),errAKR2d,'r')
legend('relerr Bernstein 2d','relerr AKR 2d')
grid
