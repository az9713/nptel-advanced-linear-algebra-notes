"""Generate the original SVG figures of the numerical laboratory page."""
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "figures"
STYLE = ("<style>text{font:15px Georgia,serif;fill:#e2e8f0}"
         ".m{fill:#94a3b8;font-size:13px}</style>")


def svg(width, height, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} '
            f'{height}" width="{width}" height="{height}" role="img">'
            f'<rect width="{width}" height="{height}" fill="#0f172a"/>'
            f'{STYLE}{body}</svg>\n')


def arrow(x1, y1, x2, y2, colour, name):
    return (f'<defs><marker id="{name}" viewBox="0 0 10 10" refX="9" '
            f'refY="5" markerWidth="7" markerHeight="7" orient="auto">'
            f'<path d="M0 0L10 5L0 10z" fill="{colour}"/></marker></defs>'
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{colour}" stroke-width="3" marker-end="url(#{name})"/>')


def projection_figure():
    """Orthogonal projection of v=(1,3) onto the line spanned by u=(3,1)."""
    u, v = (Fraction(3), Fraction(1)), (Fraction(1), Fraction(3))
    c = (v[0]*u[0] + v[1]*u[1]) / (u[0]*u[0] + u[1]*u[1])
    p = (c*u[0], c*u[1])
    r = (v[0] - p[0], v[1] - p[1])
    assert r[0]*u[0] + r[1]*u[1] == 0 and p == (Fraction(9, 5),
                                                Fraction(3, 5))
    s, ox, oy = 110, 60, 400

    def xy(point):
        return ox + float(point[0])*s, oy - float(point[1])*s
    body = ('<line x1="40" y1="400" x2="560" y2="400" stroke="#334155"/>'
            '<line x1="60" y1="30" x2="60" y2="420" stroke="#334155"/>')
    end = xy((Fraction(45, 10), Fraction(15, 10)))
    body += (f'<line x1="{ox}" y1="{oy}" x2="{end[0]:.1f}" y2="{end[1]:.1f}"'
             ' stroke="#475569" stroke-dasharray="6 5"/>')
    body += arrow(ox, oy, *xy(u), "#2dd4bf", "au")
    body += arrow(ox, oy, *xy(v), "#fb923c", "av")
    body += arrow(ox, oy, *xy(p), "#e2e8f0", "ap")
    px, py = xy(p)
    vx, vy = xy(v)
    body += (f'<line x1="{px:.1f}" y1="{py:.1f}" x2="{vx:.1f}" y2="{vy:.1f}"'
             ' stroke="#94a3b8" stroke-width="2" stroke-dasharray="4 4"/>')
    ux, uy = xy(u)
    body += (f'<text x="{ux+8:.0f}" y="{uy+18:.0f}">u = (3, 1)</text>'
             f'<text x="{vx+10:.0f}" y="{vy+4:.0f}">v = (1, 3)</text>'
             f'<text x="{px-10:.0f}" y="{py+30:.0f}">p = (9/5, 3/5)</text>'
             f'<text x="{(px+vx)/2+12:.0f}" y="{(py+vy)/2:.0f}">'
             'v − p ⟂ u</text>'
             '<text class="m" x="300" y="440">line W = span(u), dashed'
             '</text>')
    return svg(600, 460, body)


def map_figure():
    """Image of the unit square under A=[[2,1],[1,2]] with eigenlines."""
    a = ((2, 1), (1, 2))
    s, ox, oy = 70, 200, 330

    def xy(x, y):
        return ox + x*s, oy - y*s
    corners = [(0, 0), (1, 0), (1, 1), (0, 1)]
    image = [(a[0][0]*x + a[0][1]*y, a[1][0]*x + a[1][1]*y)
             for x, y in corners]
    assert image[2] == (3, 3) and image[1] == (2, 1)
    body = ('<line x1="20" y1="330" x2="580" y2="330" stroke="#334155"/>'
            '<line x1="200" y1="20" x2="200" y2="420" stroke="#334155"/>')
    for t, colour in [((1, 1), "#2dd4bf"), ((1, -1), "#fb923c")]:
        x1, y1 = xy(-1.3*t[0], -1.3*t[1])
        x2, y2 = xy(4*t[0] if t[1] > 0 else 1.3, 4*t[1] if t[1] > 0
                    else -1.3)
        body += (f'<line x1="{x1:.0f}" y1="{y1:.0f}" x2="{x2:.0f}" '
                 f'y2="{y2:.0f}" stroke="{colour}" stroke-dasharray="7 5"/>')
    square = " ".join(f"{xy(x, y)[0]:.0f},{xy(x, y)[1]:.0f}"
                      for x, y in corners)
    para = " ".join(f"{xy(x, y)[0]:.0f},{xy(x, y)[1]:.0f}" for x, y in image)
    body += (f'<polygon points="{square}" fill="#334155" '
             'fill-opacity="0.8" stroke="#94a3b8"/>'
             f'<polygon points="{para}" fill="#fb923c" fill-opacity="0.18" '
             'stroke="#fb923c" stroke-width="2"/>')
    body += ('<text x="420" y="70">eigenline x = y, λ = 3</text>'
             '<text x="300" y="420">eigenline x = −y, λ = 1</text>'
             '<text x="215" y="318" class="m">unit square</text>'
             '<text x="330" y="120">image A(square), area 3</text>')
    return svg(600, 440, body)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, text in [("projection.svg", projection_figure()),
                       ("linear-map.svg", map_figure())]:
        (OUT / name).write_text(text, encoding="utf-8", newline="\n")
    print("Wrote 2 figures")


if __name__ == "__main__":
    main()
