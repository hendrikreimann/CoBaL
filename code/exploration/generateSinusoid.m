function sinusoid = generateSinusoid(time, period, offset_time, offset_value, amplitude)
    sinusoid = sin((time - offset_time) * 1/period * 2 * pi) * amplitude + offset_value;
end