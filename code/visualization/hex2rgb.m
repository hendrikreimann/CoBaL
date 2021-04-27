function color_rgb = hex2rgb(color_hex)
    color_rgb = sscanf(color_hex(2:end),'%2x%2x%2x',[1 3])/255;
end