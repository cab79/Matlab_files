function [dp,c] = signal_detection(h,fA)

for i = 1:length(h)

    % d prime = z(h)-z(fA)
    dp(i) = norminv(h(i))-norminv(fA(i));

    % c = -0.5*[z(h)+z(fA)]
    c(i) = -0.5*(norminv(h(i))+ norminv(fA(i)));
end

