% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
% 
% 
clear all; clf;
colormap('jet');


N = 2^8;
x = linspace(-1/2,1/2,N);
[X1,X2] = meshgrid(x,x);
T = 0.1;

R0 = 0.4;

k = [0:N/2,-N/2+1:-1];
[K1,K2] = meshgrid(k,k);
Delta = -4*pi^2*(K1.^2 + K2.^2);

F_prim = @(s) (1 - 6*s).*s;
F_seconde = @(s) (1 - 12*s);

 T_vec = linspace(0,T*1.24,7)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%  1er test   %%%%%%%%%%%%%%%%%%%%%% 

epsilon = 2/N;
dist1 = abs((sqrt((X1).^2 + X2.^2 ) - R0));
dist2 = abs(X2);
dist3 = abs(X1);
dist1 = min(min(dist1,dist2),dist3);



u = 1/4*(1 - tanh(dist1/epsilon/2).^2);
R_init = sum(u(:))/N^2/(2*pi)/epsilon

 

dt =0.1*epsilon^2;
alpha =0/epsilon^2;
beta = 0/epsilon^0;
sigma = 4;
M = 1./(1 + dt*( 1*sigma*epsilon^2*Delta.^2  - Delta  +   alpha - beta*Delta));
j_sauvegarde  = 1;

 
for i=1:T/dt,
   
    
    Delta_u = (ifft2(Delta.*fft2(u)));  mu = Delta_u - F_prim(u)/epsilon^2;
    Delta_Wu = (ifft2(Delta.*fft2(F_prim(u)/epsilon^2)));
    res = sigma*epsilon^2*Delta_Wu + sigma*F_seconde(u).*(mu) + alpha*u - beta*Delta_u - F_prim(u)/epsilon^2; 
    u = real(ifft2(M.*(fft2( u + dt*res)))); 
    u = min(max(u,0),0.255);
    t1(i)=i*dt;
    Rayon1(i) =  sum(u(:))/N^2/(2*pi)/epsilon; 
    
    
    if mod(i,100)==1 
    imagesc(x,x,u);
    axis square;
    pause(0.1)
    end
    
    
     if (i*dt > T_vec(j_sauvegarde))
       
       %clf;
       imagesc(x,x,u);
       axis square
       name_title = ['t = ',num2str(i*dt)];
       title(name_title,'linewidth',2)
        
       
       name_fig = ['Test2_circle_eps1_',num2str( j_sauvegarde),'.eps'];
      
       print('-depsc', name_fig)
      
       j_sauvegarde = j_sauvegarde +1;
 

       
    end
      
end
 


