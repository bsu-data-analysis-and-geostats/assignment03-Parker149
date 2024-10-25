function ymod=gradient_search_method(A,z,v)

% equation for ice flow
ymod=80.0867-(A*(917*9.8*sin(10*pi/180)).^3).*(z.^4);  
% calc the RMSE of the function with the A paramter with raw data of v
ymod=sqrt(mean((ymod-v).^2));

