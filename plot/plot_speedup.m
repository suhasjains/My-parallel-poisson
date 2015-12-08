x = [1 2 4 8 16]
y = [1 1.99 3.93 7.73 14.5]
plot(x,y,'-o')
xlabel('Number of cores')
ylabel('Speedup')
axis([1 16 1 16]) 
title('Graph of speedup of Parallel poisson solver with Jacobi method')
