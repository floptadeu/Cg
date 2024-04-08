import tkinter as tk

# Inicializa variáveis globais
drawing_mode = 'line_bresenham'  # Modos: 'line_bresenham', 'circle_bresenham', 'line_dda'
clipping_algorithm = 'cohen_sutherland'  # Algoritmos de recorte: 'cohen_sutherland', 'liang_barsky'
points = []

def cohen_sutherland(xmin, ymin, xmax, ymax, x0, y0, x1, y1):
    """Recorta uma linha usando o algoritmo de Cohen-Sutherland."""
    INSIDE, LEFT, RIGHT, BOTTOM, TOP = 0, 1, 2, 4, 8

    def compute_outcode(x, y):
        code = INSIDE
        if x < xmin:
            code |= LEFT
        elif x > xmax:
            code |= RIGHT
        if y < ymin:
            code |= BOTTOM
        elif y > ymax:
            code |= TOP
        return code

    code0 = compute_outcode(x0, y0)
    code1 = compute_outcode(x1, y1)

    while True:
        if not (code0 | code1):
            return x0, y0, x1, y1
        elif code0 & code1:
            return None, None, None, None
        else:
            x, y = 0, 0
            code_out = code0 if code0 else code1

            if code_out & TOP:
                x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0)
                y = ymax
            elif code_out & BOTTOM:
                x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0)
                y = ymin
            elif code_out & RIGHT:
                y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0)
                x = xmax
            elif code_out & LEFT:
                y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0)
                x = xmin

            if code_out == code0:
                x0, y0 = x, y
                code0 = compute_outcode(x0, y0)
            else:
                x1, y1 = x, y
                code1 = compute_outcode(x1, y1)


def liang_barsky(xmin, ymin, xmax, ymax, x0, y0, x1, y1):
    """Recorta uma linha usando o algoritmo de Liang-Barsky."""
    dx = x1 - x0
    dy = y1 - y0
    p = [-dx, dx, -dy, dy]
    q = [x0 - xmin, xmax - x0, y0 - ymin, ymax - y0]

    t_min, t_max = 0, 1

    for i in range(4):
        if p[i] == 0:
            if q[i] < 0:
                return None, None, None, None
        else:
            t = q[i] / p[i]
            if p[i] < 0:
                t_min = max(t_min, t)
            else:
                t_max = min(t_max, t)
    
    if t_min > t_max:
        return None, None, None, None

    x0_clip = x0 + t_min * dx
    y0_clip = y0 + t_min * dy
    x1_clip = x0 + t_max * dx
    y1_clip = y0 + t_max * dy

    return x0_clip, y0_clip, x1_clip, y1_clip


def sutherland_hodgeman(subject_polygon, clip_polygon):
    """Recorta um polígono usando o algoritmo de Sutherland-Hodgeman."""
    def inside(p, p1, p2, edge):
        """Verifica se o ponto está dentro do polígono."""
        if edge == "left":
            return p[0] >= p1[0]
        elif edge == "right":
            return p[0] <= p2[0]
        elif edge == "bottom":
            return p[1] >= p1[1]
        elif edge == "top":
            return p[1] <= p2[1]

    def intersection(p1, p2, edge):
        """Calcula a interseção entre a aresta e o limite do polígono de recorte."""
        if edge == "left" or edge == "right":
            return (edge, p1[1] + (clip_x - p1[0]) * (p2[1] - p1[1]) / (p2[0] - p1[0]))
        elif edge == "bottom" or edge == "top":
            return (p1[0] + (clip_y - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1]), edge)

    output_polygon = subject_polygon

    for edge in ["left", "right", "bottom", "top"]:
        clip_x, clip_y = clip_polygon[0]
        input_polygon = output_polygon
        output_polygon = []

        for i in range(len(input_polygon)):
            p1 = input_polygon[i]
            p2 = input_polygon[(i + 1) % len(input_polygon)]

            if inside(p1, p1, p2, edge) and inside(p2, p1, p2, edge):
                output_polygon.append(p2)
            elif inside(p1, p1, p2, edge) and not inside(p2, p1, p2, edge):
                output_polygon.append(intersection(p1, p2, edge))
            elif not inside(p1, p1, p2, edge) and inside(p2, p1, p2, edge):
                output_polygon.append(intersection(p1, p2, edge))
                output_polygon.append(p2)

    return output_polygon


def set_pixel(canvas, x, y, color='red'):
    """Desenha um pixel no canvas."""
    canvas.create_rectangle(x, y, x+1, y+1, outline=color, fill=color)

def line_dda(canvas, x0, y0, xEnd, yEnd):
    """Desenha uma linha usando o algoritmo DDA."""
    dx = xEnd - x0
    dy = yEnd - y0
    steps = max(abs(dx), abs(dy))
    x_increment = dx / steps
    y_increment = dy / steps
    x = x0
    y = y0
    
    for _ in range(int(steps)):
        set_pixel(canvas, round(x), round(y))
        x += x_increment
        y += y_increment

def line_bresenham(canvas, x0, y0, xEnd, yEnd):
    """Desenha uma linha usando o algoritmo de Bresenham."""
    dx = abs(xEnd - x0)
    dy = abs(yEnd - y0)
    sx = 1 if xEnd - x0 > 0 else -1
    sy = 1 if yEnd - y0 > 0 else -1
    err = dx - dy

    while True:
        set_pixel(canvas, x0, y0)
        if x0 == xEnd and y0 == yEnd:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x0 += sx
        if e2 < dx:
            err += dx
            y0 += sy

def boundary_fill(canvas, x, y, boundary_color, fill_color):
    current_color = canvas.get_pixel(x, y)
    if current_color != boundary_color and current_color != fill_color:
        canvas.set_pixel(x, y, fill_color)
        boundary_fill(canvas, x + 1, y, boundary_color, fill_color)
        boundary_fill(canvas, x - 1, y, boundary_color, fill_color)
        boundary_fill(canvas, x, y + 1, boundary_color, fill_color)
        boundary_fill(canvas, x, y - 1, boundary_color, fill_color)
def flood_fill(canvas, x, y, target_color, fill_color):
    current_color = canvas.get_pixel(x, y)
    if current_color == target_color:
        canvas.set_pixel(x, y, fill_color)
        flood_fill(canvas, x + 1, y, target_color, fill_color)
        flood_fill(canvas, x - 1, y, target_color, fill_color)
        flood_fill(canvas, x, y + 1, target_color, fill_color)
        flood_fill(canvas, x, y - 1, target_color, fill_color)
def scanline_fill(canvas, edges, fill_color):
    ymin = min(edges, key=lambda edge: edge[0])[0]
    ymax = max(edges, key=lambda edge: edge[1])[1]

    for y in range(ymin, ymax + 1):
        intersections = []
        for edge in edges:
            if edge[0] <= y < edge[1]:
                intersections.append(edge[2] + (y - edge[0]) / (edge[1] - edge[0]) * (edge[3] - edge[2]))

        intersections.sort()
        for i in range(0, len(intersections), 2):
            for x in range(int(intersections[i]), int(intersections[i + 1]) + 1):
                canvas.set_pixel(x, y, fill_color)


def plot_point(canvas, xc, yc, x, y):
    """Desenha pontos em todas as oitavas a partir de um centro."""
    points = [
        (xc+x, yc+y), (xc+x, yc-y), (xc+y, yc+x), (xc+y, yc-x),
        (xc-x, yc-y), (xc-y, yc-x), (xc-x, yc+y), (xc-y, yc+x)
    ]
    for point in points:
        set_pixel(canvas, *point)

def circle_bresenham(canvas, xc, yc, r):
    """Desenha um círculo usando o algoritmo de Bresenham."""
    x, y = 0, r
    pk = 5/4 - r
    plot_point(canvas, xc, yc, x, y)
    while x < y:
        x += 1
        if pk < 0:
            pk += 2*x + 1
        else:
            y -= 1
            pk += 2*(x - y) + 1
        plot_point(canvas, xc, yc, x, y)

def on_click(event):
    """Manipula cliques no canvas."""
    global points
    canvas = event.widget
    point = (event.x, event.y)
    points.append(point)
    if drawing_mode == 'circle_bresenham' and len(points) == 1:
        # Para círculos, o primeiro clique define o centro, o segundo define o raio
        return
    elif len(points) == 2:
        if drawing_mode == 'line_dda':
            x0, y0, x1, y1 = points[0][0], points[0][1], points[1][0], points[1][1]
            if clipping_algorithm == 'liang_barsky':
                x0_clip, y0_clip, x1_clip, y1_clip = liang_barsky(0, 0, canvas.winfo_width(), canvas.winfo_height(), x0, y0, x1, y1)
            else:
                x0_clip, y0_clip, x1_clip, y1_clip = cohen_sutherland(0, 0, canvas.winfo_width(), canvas.winfo_height(), x0, y0, x1, y1)
            if x0_clip is not None and y0_clip is not None and x1_clip is not None and y1_clip is not None:
                line_dda(canvas, x0_clip, y0_clip, x1_clip, y1_clip)
        elif drawing_mode == 'line_bresenham':
            x0, y0, x1, y1 = points[0][0], points[0][1], points[1][0], points[1][1]
            if clipping_algorithm == 'liang_barsky':
                x0_clip, y0_clip, x1_clip, y1_clip = liang_barsky(0, 0, canvas.winfo_width(), canvas.winfo_height(), x0, y0, x1, y1)
            else:
                x0_clip, y0_clip, x1_clip, y1_clip = cohen_sutherland(0, 0, canvas.winfo_width(), canvas.winfo_height(), x0, y0, x1, y1)
            if x0_clip is not None and y0_clip is not None and x1_clip is not None and y1_clip is not None:
                line_bresenham(canvas, x0_clip, y0_clip, x1_clip, y1_clip)
        elif drawing_mode == 'circle_bresenham':
            r = int(((points[0][0] - points[1][0]) ** 2 + (points[0][1] - points[1][1]) ** 2) ** 0.5)
            circle_bresenham(canvas, *points[0], r)
        points.clear()

def set_drawing_mode(mode):
    """Define o modo de desenho."""
    global drawing_mode
    drawing_mode = mode
    points.clear()

def set_clipping_algorithm(algorithm):
    """Define o algoritmo de recorte."""
    global clipping_algorithm
    clipping_algorithm = algorithm

def clear_canvas():
    """Limpa o canvas."""
    canvas.delete("all")
    points.clear()

def zoom(event):
    """Função para realizar o zoom usando o scroll do mouse."""
    global scale_factor
    if event.delta > 0:
        canvas.scale("all", event.x, event.y, 1.1, 1.1)
    elif event.delta < 0:
        canvas.scale("all", event.x, event.y, 0.9, 0.9)

def move_canvas(event):
    """Função para mover o canvas usando o botão contrário do mouse."""
    canvas.scan_mark(event.x, event.y)
    canvas.bind("<B2-Motion>", lambda e: canvas.scan_dragto(e.x, e.y, gain=1))

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Desenhando com Algoritmos")

    # Criando o canvas
    canvas = tk.Canvas(root, width=600, height=400, bg="white")
    canvas.pack(expand=True, fill="both")
    canvas.bind("<Button-1>", on_click)
    canvas.bind("<MouseWheel>", zoom)
    canvas.bind("<Button-3>", move_canvas)

    # Botões para selecionar o modo de desenho
    line_bresenham_btn = tk.Button(root, text="Linha (Bresenham)", command=lambda: set_drawing_mode('line_bresenham'))
    line_bresenham_btn.pack(side=tk.LEFT, padx=5, pady=5)
    circle_bresenham_btn = tk.Button(root, text="Círculo (Bresenham)", command=lambda: set_drawing_mode('circle_bresenham'))
    circle_bresenham_btn.pack(side=tk.LEFT, padx=5, pady=5)
    line_dda_btn = tk.Button(root, text="Linha (DDA)", command=lambda: set_drawing_mode('line_dda'))
    line_dda_btn.pack(side=tk.LEFT, padx=5, pady=5)

    # Botões para selecionar o algoritmo de recorte
    cohen_sutherland_btn = tk.Button(root, text="Cohen-Sutherland", command=lambda: set_clipping_algorithm('cohen_sutherland'))
    cohen_sutherland_btn.pack(side=tk.LEFT, padx=5, pady=5)
    liang_barsky_btn = tk.Button(root, text="Liang-Barsky", command=lambda: set_clipping_algorithm('liang_barsky'))
    liang_barsky_btn.pack(side=tk.LEFT, padx=5, pady=5)

    # Botão para limpar o canvas
    clear_btn = tk.Button(root, text="Limpar", command=clear_canvas)
    clear_btn.pack(side=tk.LEFT, padx=5, pady=5)

    root.mainloop()
