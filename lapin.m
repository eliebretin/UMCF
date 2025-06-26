function donne = lapin()
fiche=fopen('lapin.dat','r') ;
N = 128;
donne = zeros(N,N,N);

for i=1:128,
donne(:,:,i)=fscanf(fiche,'%g %g',[N,N]);
end

fclose(fiche);

donne = permute(donne,[3 2 1]);

% x=linspace(0,1,N);
% y=linspace(0,1,N);
% z=linspace(0,1,N);
% v=donne;
% 
%  p = patch(isosurface(x,y,z,v,0));
%  isonormals(x,y,z,v,p)
%  set(p,'FaceColor','cyan','EdgeColor','none');
%  daspect([1 1 1])
%  view(3); axis tight
%  camlight
%  lighting gouraud
%  pause(0.1)



end