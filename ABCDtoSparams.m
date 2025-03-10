function [S11, S12, S21, S22] = ABCDtoSparams(A, B, C, D, Z0)

S11 = (A + B./Z0 - C.*Z0 - D)./ ...
    (A + B./Z0 + C.*Z0 + D);
S12 = 2.*(A.*D - B.*C)./ ...
    (A + B./Z0 + C.*Z0 + D);
S21 = 2./(A + B./Z0 + C.*Z0 + D);
S22 = (-A + B./Z0 - C.*Z0 + D)./ ...
    (A + B./Z0 + C.*Z0 + D);

end