function [vec_rate] = compute_exp_rate_conv(x_bar, n_iter)
vec_rate = zeros(1,n_iter);
for i=3:(length(x_bar)-1)
    e_succ = norm(x_bar(i+1,:)-x_bar(i,:));
    disp("xbar")
    disp(x_bar(i-2,:))
    disp(x_bar(i-1,:))
    disp(x_bar(i,:))
    disp(x_bar(i+1,:))
    e_k = norm(x_bar(i,:)-x_bar(i-1,:));
    e_prev = norm(x_bar(i-1,:)-x_bar(i-2,:));
    % excluded iterations with the same x_bar in order to avoid NAN
    % solutions
    if e_k ~= 0 && e_prev ~= 0 && e_succ ~= 0
        num = log(e_succ/e_k);
        den = log(e_k/e_prev);
        % excluded iterations with the same increments (e_k = e_prev)
        % solutions
        if den ~= 0 || num*den > 0  %perchè escono fuori cose non maggiori di 0?
            vec_rate(i) = num/den;
        end
    end
end
vec_rate = vec_rate(3:(length(x_bar)-1));
plot(1:length(vec_rate), vec_rate)