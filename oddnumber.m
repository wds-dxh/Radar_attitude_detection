function n = oddnumber(n)
    % 确保输入为奇数，如果不是，则减去1使之成为奇数
    if mod(n, 2) == 0
        n = n - 1;
    end
end