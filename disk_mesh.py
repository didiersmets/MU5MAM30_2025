from manim import *
import numpy as np

class DiskMesh(Scene):
    def construct(self):
        # Parameters for the mesh
        N = 4
        R = 3.0
        
        # 1. Build Vertices
        vertices = []
        vertices.append(np.array([0.0, 0.0, 0.0]))
        
        for i in range(1, N + 1):
            for j in range(6 * i):
                angle = PI * j / (3.0 * i)
                r = R * i / N
                vertices.append(np.array([r * np.cos(angle), r * np.sin(angle), 0.0]))
                
        # 2. Build Triangles (using exact C++ indexing logic)
        triangles = []
        
        def set_tri(a, b, c):
            triangles.append((a, b, c))

        # Initial six triangles around center
        set_tri(1, 2, 0)
        set_tri(2, 3, 0)
        set_tri(3, 4, 0)
        set_tri(4, 5, 0)
        set_tri(5, 6, 0)
        set_tri(6, 1, 0)

        # Further rings
        for i in range(2, N + 1):
            aux = 0
            for j in range(6 * i - 1):
                if j % i != 0:
                    set_tri(2 + 3 * (i - 2) * (i - 1) + aux,
                            1 + 3 * (i - 1) * (i - 2) + aux,
                            1 + 3 * (i - 1) * i + j)
                    set_tri(2 + 3 * (i - 2) * (i - 1) + aux,
                            1 + 3 * (i - 1) * i + j,
                            2 + 3 * (i - 1) * i + j)
                    aux += 1
                else:
                    set_tri(1 + 3 * (i - 1) * i + j,
                            2 + 3 * (i - 1) * i + j,
                            1 + 3 * (i - 2) * (i - 1) + aux)
                            
            set_tri(1 + 3 * (i - 1) * i,
                    3 * (i - 1) * i + 6 * i,
                    2 + aux + (i - 1) * (3 * (i - 2) - 6))
            set_tri(2 + aux + (i - 1) * (3 * (i - 2) - 6),
                    1 + 3 * (i - 2) * (i - 1) + aux,
                    3 * (i - 1) * i + 6 * i)

        # 3. Extract unique edges from triangles for wireframe rendering
        unique_edges = set()
        for t in triangles:
            unique_edges.add(tuple(sorted([t[0], t[1]])))
            unique_edges.add(tuple(sorted([t[1], t[2]])))
            unique_edges.add(tuple(sorted([t[2], t[0]])))

        # 4. Create Manim VGroups
        dots = VGroup()
        lines = VGroup()
        
        # Calculate boundary start index based on your C++ logic
        vtx_count = len(vertices)
        boundary_start = 1 + 3 * (N - 1) * N

        # Create vertices (Dots)
        for idx, pos in enumerate(vertices):
            # Highlight boundary nodes in RED, inner nodes in BLUE
            color = RED if idx >= boundary_start else BLUE
            dot = Dot(point=pos, radius=0.06, color=color)
            dots.add(dot)

        # Create edges (Lines)
        for edge in unique_edges:
            start_pos = vertices[edge[0]]
            end_pos = vertices[edge[1]]
            line = Line(start_pos, end_pos, stroke_width=2, color=WHITE).set_opacity(0.5)
            lines.add(line)

        # 5. Animation sequence
        self.play(Write(Text("Generating Disk Mesh").to_edge(UP)))
        
        # Animate inner vertices and rings growing outward
        self.play(Create(dots), run_time=3)
        self.wait(0.5)
        
        # Draw the wireframe connections
        self.play(Create(lines), run_time=4)
        self.wait(2)