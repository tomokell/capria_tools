function outcol = tofadecolor(col,p)
% Fade a 3-component color so it's closer to white by a proportion p
    w = [1 1 1];
    outcol = (w - col(:)')*p + col(:)';
    outcol =reshape(outcol,size(col));
end

