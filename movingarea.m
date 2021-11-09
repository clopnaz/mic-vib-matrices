function t = movingarea(g)
    t(1,1) = sym('0');
    t(1,2) = 0;
    t(2,1) = 0;
    t(2,2) = 1-g;
end