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
dist1 = (sqrt((X1).^2 + X2.^2 ) - R0);
u = 1/4*(1 - tanh(dist1/epsilon/2).^2);
R_init = sum(u(:))/N^2/(2*pi)/epsilon

 

dt =0.01*epsilon^2;
alpha =0/epsilon^2;
beta = 0/epsilon^0;
sigma = 2;
M = 1./(1 + dt*( 1*sigma*epsilon^2*Delta.^2  - Delta  +   alpha - beta*Delta));
j_sauvegarde  = 1;

 
for i=1:T/dt,
   
    
    Delta_u = (ifft2(Delta.*fft2(u)));  mu = Delta_u - F_prim(u)/epsilon^2;
    Delta_Wu = (ifft2(Delta.*fft2(F_prim(u)/epsilon^2)));
    res = sigma*epsilon^2*Delta_Wu + sigma*F_seconde(u).*(mu) + alpha*u - beta*Delta_u - F_prim(u)/epsilon^2; 
    u = real(ifft2(M.*(fft2( u + dt*res)))); 
    u = min(max(u,0),0.25);
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
        
       
       name_fig = ['Test_circle_eps1_',num2str( j_sauvegarde),'.eps'];
      
       print('-depsc', name_fig)
      
       j_sauvegarde = j_sauvegarde +1;
 

       
    end
      
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù

epsilon = 3/N;
dist1 = (sqrt((X1).^2 + X2.^2 ) - R0);
u = 1/4*(1 - tanh(dist1/epsilon/2).^2);
R_init = sum(u(:))/N^2/(2*pi)/epsilon

 

dt =0.01*epsilon^2;
T = 0.1;
alpha =0/epsilon^2;
beta = 0/epsilon^0;
sigma = 1;
M = 1./(1 + dt*( 1*sigma*epsilon^2*Delta.^2  - Delta  +   alpha - beta*Delta));
j_sauvegarde  = 1;

 
for i=1:T/dt,
   
    
    Delta_u = (ifft2(Delta.*fft2(u)));  mu = Delta_u - F_prim(u)/epsilon^2;
    Delta_Wu = (ifft2(Delta.*fft2(F_prim(u)/epsilon^2)));
    res = sigma*epsilon^2*Delta_Wu + sigma*F_seconde(u).*(mu) + alpha*u - beta*Delta_u - F_prim(u)/epsilon^2; 
    u = real(ifft2(M.*(fft2( u + dt*res)))); 
    u = min(max(u,0),0.25);
    t2(i)=i*dt;
    Rayon2(i) =  sum(u(:))/N^2/(2*pi)/epsilon; 
    
    
    if mod(i,100)==1 
    imagesc(u);
    axis square;
    pause(0.1)
    end
    
    
     if (i*dt > T_vec(j_sauvegarde))
       
      % clf;
       imagesc(x,x,u);
         axis square
       name_title = ['t = ',num2str(i*dt),'\epsilon = 3/N'];
       title(name_title,'linewidth',2)

       
       name_fig = ['Test_circle_eps2_2',num2str( j_sauvegarde),'.eps'];
      
       print('-depsc', name_fig)
      
       j_sauvegarde = j_sauvegarde +1;
 

       
    end
      
end








clf;

%%%%%%%%%%%%% Resolution du flot %%%%%%%%%%%%
f = @(R) -1/R + 1/2*epsilon^3*1/R^3*sigma;
R = R0;
for i=1:T/dt;
   k = R + dt*f(R);
   R = max(R + dt*(f(R)+f(k))/2,0);
   Rayon_th(i)=R;
    
end


plot(t1,sqrt(max(R0^2 - 2*t1,0)),'r','LineWidth', 2);
hold on;
plot(t1,Rayon1,'g','LineWidth', 2);
plot(t2,Rayon2,'b','LineWidth', 2);
legend('Exact','\epsilon = 1.5/N','\epsilon = 3/N')

name_title = ['Radius R(t) for different values of \epsilon'];
title(name_title,'linewidth',2)

name_fig = ['Test_circle__erreur_eps.eps'];
print('-depsc', name_fig)
% 
 