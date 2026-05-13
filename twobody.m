function ydot = twobody(t, y, mu)
    r = y(1:3);
    v = y(4:6);
    ydot = [v; -mu/norm(r)^3*r];
end