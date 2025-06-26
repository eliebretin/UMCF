% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
% 
% 
N = 2^8;
epsilon = 3/N;
x = linspace(-1/2,1/2,N);
[X1,X2] = meshgrid(x,x);
R0 = 0.4;
dist1 = (sqrt((X1).^2 + X2.^2 ) - R0);
u = 1/4*(1 - tanh(dist1/epsilon/2).^2);

R_init = sum(u(:))/N^2/(2*pi)/epsilon

 k = [0:N/2,-N/2+1:-1];
 [K1,K2] = meshgrid(k,k);
 Delta = -4*pi^2*(K1.^2 + K2.^2);
 dt =1*epsilon^2;
 T = 10;
 
 
 W = @(s) (1 - 6*s).*s;
 W_prim = @(s) (1 - 12*s);
 
 
 
 alpha =2/epsilon^2;
 beta = 2/epsilon^0;
 sigma = 2;

 M = 1./(1 + dt*( 1*sigma*epsilon^2*Delta.^2  - Delta  +   alpha - beta*Delta));

 
 colormap('jet');
 
 u = min(max(u,0),0.25);
 
 
 
  
for i=1:T/dt,
    
    
    Delta_u = (ifft2(Delta.*fft2(u)));  mu = Delta_u - W(u)/epsilon^2;
    Delta_Wu = (ifft2(Delta.*fft2(W(u)/epsilon^2)));
    res = sigma*epsilon^2*Delta_Wu + sigma*W_prim(u).*(mu) + alpha*u - beta*Delta_u - W(u)/epsilon^2; 
    u = real(ifft2(M.*(fft2( u + dt*res)))); 
    u = min(max(u,0),0.25);
    
    
    
    
    if mod(i,100)==1 
    imagesc(u);
    axis square;
    pause(0.1)
    end
      
 end
% 
 