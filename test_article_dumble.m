% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
% 
% 
clear all; clf;
colormap('jet');


N = 2^7;
x = linspace(-1/2,1/2,N);
[X1,X2,X3] = meshgrid(x,x,x);

R0 = 0.3;


k = [0:N/2,-N/2+1:-1];
[K1,K2,K3] = meshgrid(k,k,k);
Delta = -4*pi^2*(K1.^2 + K2.^2 + K3.^2);

F_prim = @(s) (1 - 6*s).*s;
F_seconde = @(s) (1 - 12*s);

T = 0.012;
T_vec = linspace(0,T*1.01,10)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%  1er test   %%%%%%%%%%%%%%%%%%%%%% 

epsilon = 2/N;
R_1 = 0.10; R_2 = 0.3; R_3 = 0.15;
dist = max((sqrt(X1.^2 +X2.^2 )-R_1), X3.^2 - R_2^2);
dist = min(min(dist,sqrt(0.5*X1.^2 + 0.5*X2.^2 + (X3-R_2).^2)-R_3),sqrt(0.5*X1.^2 + 0.5*X2.^2 + (X3+R_2).^2)-R_3);
u = 1/4*(1 - tanh(dist/epsilon/2).^2);
R_init = sum(u(:))/N^2/(2*pi)/epsilon

affiche_solution_3d2(x,4*u,0*u);
view(-100,20);
axis([-0.5,0.5,-0.5,0.5,-0.5,0.5])


dt =0.01*epsilon^2;
alpha =0/epsilon^2;
beta = 0/epsilon^0;
sigma = 2;
M = 1./(1 + dt*( 1*sigma*epsilon^2*Delta.^2  - Delta  +   alpha - beta*Delta));
j_sauvegarde  = 1;

 
for i=1:T/dt,
    i
    
    Delta_u = (ifftn(Delta.*fftn(u)));  mu = Delta_u - F_prim(u)/epsilon^2;
    Delta_Wu = (ifftn(Delta.*fftn(F_prim(u)/epsilon^2)));
    res = sigma*epsilon^2*Delta_Wu + sigma*F_seconde(u).*(mu) + alpha*u - beta*Delta_u - F_prim(u)/epsilon^2; 
    u = real(ifftn(M.*(fftn( u + dt*res)))); 
    u = min(max(u,0),0.25);
  
   
    
    if mod(i,100)==1 
     clf;
       
       affiche_solution_3d2(x,3*u,0*u);
      view(-100,20);
      axis([-0.5,0.5,-0.5,0.5,-0.5,0.5])
      pause(0.1)
    end
    
    
     if (i*dt > T_vec(j_sauvegarde))
       
       clf;
       
       affiche_solution_3d2(x,3*u,0*u);
      view(-100,20);
      axis([-0.5,0.5,-0.5,0.5,-0.5,0.5])

       
       axis square
       name_title = ['t = ',num2str(i*dt),'\epsilon = 2/N'];
       title(name_title,'linewidth',2)
        
       
       name_fig = ['Test1_dumble_',num2str( j_sauvegarde),'.eps'];
      
       print('-depsc', name_fig)
      
       j_sauvegarde = j_sauvegarde +1;
 
       pause(0.1)
       
    end
      
end






 
%%%%%%%%%%%%%%%%%%%%%%%%%% 2 test   %%%%%%%%%%%%%%%%%%%%%% 

epsilon = 1/N;
R_1 = 0.10; R_2 = 0.3; R_3 = 0.15;
dist = max((sqrt(X1.^2 +X2.^2 )-R_1), X3.^2 - R_2^2);
dist = min(min(dist,sqrt(0.5*X1.^2 + 0.5*X2.^2 + (X3-R_2).^2)-R_3),sqrt(0.5*X1.^2 + 0.5*X2.^2 + (X3+R_2).^2)-R_3);
u = 1/4*(1 - tanh(dist/epsilon/2).^2);
R_init = sum(u(:))/N^2/(2*pi)/epsilon

affiche_solution_3d2(x,4*u,0*u);
view(-100,20);
axis([-0.5,0.5,-0.5,0.5,-0.5,0.5])


dt =0.01*epsilon^2;
alpha =0/epsilon^2;
beta = 0/epsilon^0;
sigma = 2;
M = 1./(1 + dt*( 1*sigma*epsilon^2*Delta.^2  - Delta  +   alpha - beta*Delta));
j_sauvegarde  = 1;

 
for i=1:T/dt,
    i
    
    Delta_u = (ifftn(Delta.*fftn(u)));  mu = Delta_u - F_prim(u)/epsilon^2;
    Delta_Wu = (ifftn(Delta.*fftn(F_prim(u)/epsilon^2)));
    res = sigma*epsilon^2*Delta_Wu + sigma*F_seconde(u).*(mu) + alpha*u - beta*Delta_u - F_prim(u)/epsilon^2; 
    u = real(ifftn(M.*(fftn( u + dt*res)))); 
    u = min(max(u,0),0.25);
  
   
    
    if mod(i,100)==1 
     clf;
       
       affiche_solution_3d2(x,3*u,0*u);
      view(-100,20);
      axis([-0.5,0.5,-0.5,0.5,-0.5,0.5])
      pause(0.1)
    end
    
    
     if (i*dt > T_vec(j_sauvegarde))
       
       clf;
       
       affiche_solution_3d2(x,3*u,0*u);
      view(-100,20);
      axis([-0.5,0.5,-0.5,0.5,-0.5,0.5])

       
       axis square
       name_title = ['t = ',num2str(i*dt),'\epsilon = 1/N'];
       title(name_title,'linewidth',2)
        
       
       name_fig = ['Test2_dumble_',num2str( j_sauvegarde),'.eps'];
      
       print('-depsc', name_fig)
      
       j_sauvegarde = j_sauvegarde +1;
 
       pause(0.1)
       
    end
      
end
 






 


