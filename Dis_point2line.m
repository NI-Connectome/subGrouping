function d = dis_point2line(x,y,k,c)

%kx-y+c=0

d = abs((k*x-y+c)/sqrt(k^2+1));