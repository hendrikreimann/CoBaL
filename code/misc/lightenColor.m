function lighterColor = lightenColor(color, ratio)
    lighterColor = color + (ones(size(color)) - color) * ratio;

end