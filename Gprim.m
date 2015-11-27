function [out2] = Gprim(arg2,a2)

    out2 = a2 - a2 * tanh(a2 * arg2).^2;

end