function [vec_rate] = compute_exp_rate_conv2(x_bar, n_iter, x_opt)
vec_rate = zeros(1,n_iter);
for i=1:(length(x_bar)-1)
    e_succ = norm(x_bar(i+1,:)-x_opt);
    e_k = norm(x_bar(i,:)-x_opt);
    if e_k ~= 0
        ratio = e_succ/e_k;
        vec_rate(i) = ratio;
    end
end
vec_rate = vec_rate(1:(length(x_bar)-1));
plot(1:length(vec_rate), vec_rate)