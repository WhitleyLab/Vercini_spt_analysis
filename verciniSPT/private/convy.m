function y = convy
    p = get(gcf,'currentPoint');
    fs = get(gcf,'position');
    as = get(gca,'position');
    aywp = as(4)*fs(4);
    ayyo = as(2)*fs(4);
    ayp = p(2)-ayyo;
    yl = get(gca,'ylim');
    y = ayp/aywp * (yl(2)-yl(1)) + yl(1);
end