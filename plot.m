data = load('data');
k=1;
for (i=1:102)
  for (j=1:102)
    z(i,j)=data(k,3);
    k=k+1;
  end
end
surf (z)
  