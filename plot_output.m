data = load('data0');
k=1;

z = reshape(data(:,3),1000,1000);

surf (z)

%mesh (data(:,1),data(:,2),data(:,3))
 
%mesh(data(:,3)) 
