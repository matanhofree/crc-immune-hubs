function v = hexColor(str)
    r = hex2dec(str(1:2))/255;
    g = hex2dec(str(3:4))/255;
    b = hex2dec(str(5:6))/255;
    v = [r,g,b]';
end