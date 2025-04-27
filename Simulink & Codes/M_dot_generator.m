function y = M_dot_generator(M,theta,theta_dot)
[m n] = size(M);

for i = 1:m
    for j = 1:n
      y(i,j) =   jacobian(M(i,j),theta) * theta_dot;
    end
end

y = simplify(y);