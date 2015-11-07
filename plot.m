data = load('data5');
k=1;

z = reshape(data(:,3),50,50);

surf (z)

%mesh (data(:,1),data(:,2),data(:,3))
 
%mesh(data(:,3)) 
