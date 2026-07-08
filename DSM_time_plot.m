dim_=[0 6 12 24 36 48 60 72 96];
dim_time=[0 8.410 16.39 41.06 88.64 168.89 277.05 432.26 912.80]
dim_i =0:100
plot(dim_,dim_time,'o-')
grid on 
hold on
plot(dim_i,0.127*dim_i.^2-3.098*dim_i+22.52,'r-')
