function x = convx
    p = get(gcf,'currentPoint');
    fs = get(gcf,'position');
    as = get(gca,'position');
    axwp = as(3)*fs(3);
    axxo = as(1)*fs(3);
    axp = p(1)-axxo;
    xl = get(gca,'xlim');
    x = axp/axwp * (xl(2)-xl(1)) + xl(1);
end