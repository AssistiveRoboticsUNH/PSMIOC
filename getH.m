function H = getH(tmp, W)
    H = zeros(size(tmp{1}));
    for i = 1:length(W)
        H = H + W(i)*tmp{i};
    end
end

