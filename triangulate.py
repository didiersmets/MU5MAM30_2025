from manim import *

class TriangulationDemo(MovingCameraScene):
    def construct(self):
        # --- Problem Parameters ---
        self.N = 3  # Number of vertices per side (N*N grid)
        self.a = 5.0 # Width
        self.b = 5.0 # Height
        self.corner_coords = np.array([[-self.a/2, -self.b/2, 0], [self.a/2, self.b/2, 0]])

        # Create title
        self.title = Text("Triangulation & Boundary Marking", font_size=36).to_edge(UP)
        self.add(self.title)

        # --- STEP 1: Generate Grid and Indices ---
        self.generate_grid()
        self.add_vertex_labels()
        self.wait(2)

        # --- STEP 2: Zoom and Show Triangle Order ---
        # Focus on the bottom-left square: vertices (0, 1, N, N+1) -> (0, 1, 3, 4)
        self.zoom_in_to_cell(target_vertex=0)
        self.show_winding_order()
        self.wait(3)
        self.zoom_out()

        # --- STEP 3: Complete Grid Triangulation ---
        self.triangulate_rest()
        self.wait(2)

        # --- STEP 4: Boundary Vertex Marking ---
        self.mark_corners()
        self.wait(2)
        self.mark_edges()
        self.wait(3)

    # --- HELPER FUNCTIONS ---

    def generate_grid(self):
        # Dictionaries to store dots and label objects
        self.dots = {}
        self.vertex_group = VGroup()

        # Iterate as N*N (row by row, left to right, bottom-to-top coordinates)
        # i is row (y), j is column (x)
        for i in range(self.N):
            for j in range(self.N):
                idx = self.N * i + j
                
                # Formula: x = a*(j / (N-1)), y = b*(i / (N-1)), center the grid
                x = (self.a * j / (self.N - 1)) - (self.a / 2)
                y = (self.b * i / (self.N - 1)) - (self.b / 2)
                pos = np.array([x, y, 0])
                
                dot = Dot(pos, radius=0.1, color=WHITE)
                self.dots[idx] = dot
                self.vertex_group.add(dot)

        step1_text = Text("1. Grid of 3x3 Vertices Generated", font_size=24, color=BLUE).next_to(self.title, DOWN)
        self.play(Write(step1_text), FadeIn(self.vertex_group))
        self.wait(1)
        self.remove(step1_text)

    def add_vertex_labels(self):
        self.labels = {}
        self.label_group = VGroup()
        for idx, dot in self.dots.items():
            # Label row-by-row logic (bottom-to-top, left-to-right)
            label = Text(str(idx), font_size=20).next_to(dot, DR, buff=0.1)
            self.labels[idx] = label
            self.label_group.add(label)
        
        step2_text = Text("2. Row-by-Row Vertex Indices Added", font_size=24, color=YELLOW).next_to(self.title, DOWN)
        self.play(Write(step2_text))
        self.play(FadeIn(self.label_group), run_time=1.5)
        self.wait(1)
        self.remove(step2_text)

    def zoom_in_to_cell(self, target_vertex):
        step3_text = Text("3. Zooming into Cell for Triangulation", font_size=24, color=WHITE).next_to(self.title, DOWN)
        self.add(step3_text)

        # Scale the whole camera and shift
        focus_center = self.dots[target_vertex].get_center() + np.array([0.5, 0.5, 0]) # Move to center of cell
        self.play(
            self.camera.frame.animate.scale(0.3).move_to(focus_center),
            # Fade titles out for cleaner detail view
            FadeOut(self.title),
            FadeOut(step3_text),
            run_time=2
        )

    def show_winding_order(self):
        # Indices of bottom-left cell: (0, 1, 3, 4)
        v0, v1, v2, v3 = 0, 1, 3, 4 
        
        # Triangle 0 (bottom-right): (v, v+1, v+1+N) -> (0, 1, 4)
        winding0_indices = [v0, v1, v3]
        t0_label = Text("Tri 0:\n(0, 1, 4)", font_size=16, color=BLUE_C).next_to(self.dots[v1], LEFT, buff=0.4).shift(UP*0.2)
        self.play(Write(t0_label))
        
        # Explicit lines, specific vertex indices asked for
        tri0_poly = Polygon(self.dots[v0].get_center(), self.dots[v1].get_center(), self.dots[v3].get_center(), color=BLUE_C, fill_opacity=0.4, stroke_width=2)
        tri0_labels = VGroup(
            Text("0", font_size=16, color=BLUE_C).next_to(self.dots[v0], DL, buff=0.1),
            Text("1", font_size=16, color=BLUE_C).next_to(self.dots[v1], DR, buff=0.1),
            Text("2", font_size=16, color=BLUE_C).next_to(self.dots[v3], UR, buff=0.1)
        )
        self.play(Create(tri0_poly))
        self.play(FadeIn(tri0_labels))
        self.wait(1)

        # Triangle 1 (top-left): (v, v+1+N, v+N) -> (0, 4, 3)
        winding1_indices = [v0, v3, v2]
        t1_label = Text("Tri 1:\n(0, 4, 3)", font_size=16, color=GREEN_C).next_to(self.dots[v3], LEFT, buff=0.4).shift(UP*0.2)
        self.play(Write(t1_label))

        tri1_poly = Polygon(self.dots[v0].get_center(), self.dots[v3].get_center(), self.dots[v2].get_center(), color=GREEN_C, fill_opacity=0.4, stroke_width=2)
        tri1_labels = VGroup(
            Text("3", font_size=16, color=GREEN_C).next_to(self.dots[v0], UL, buff=0.1),
            Text("4", font_size=16, color=GREEN_C).next_to(self.dots[v3], UR, buff=0.1),
            Text("5", font_size=16, color=GREEN_C).next_to(self.dots[v2], UL, buff=0.1)
        )
        self.play(Create(tri1_poly))
        self.play(FadeIn(tri1_labels))
        self.wait(1)

        # Store for cleanup
        self.cell_visuals = VGroup(t0_label, t1_label, tri0_poly, tri1_poly, tri0_labels, tri1_labels)

    def zoom_out(self):
        self.play(
            self.camera.frame.animate.scale(1/0.3).move_to(ORIGIN),
            FadeOut(self.cell_visuals),
            FadeIn(self.title), # Restore title
            run_time=1.5
        )

    def triangulate_rest(self):
        # Title of step
        step4_text = Text("4. Applying Logic to Remaining Cells", font_size=24, color=WHITE).next_to(self.title, DOWN)
        self.play(Write(step4_text))

        self.all_tris = VGroup()
        for i in range(self.N - 1): # i is row index of cell
            for j in range(self.N - 1): # j is column index of cell
                # row major index of the bottom-left vertex of the current cell
                v = self.N * i + j 
                
                # Skip the bottom-left square we already detailed visually
                if v == 0: continue

                t0 = Polygon(self.dots[v].get_center(), self.dots[v+1].get_center(), self.dots[v+1+self.N].get_center(), color=BLUE_E, fill_opacity=0.2, stroke_width=2)
                t1 = Polygon(self.dots[v].get_center(), self.dots[v+1+self.N].get_center(), self.dots[v+self.N].get_center(), color=GREEN_E, fill_opacity=0.2, stroke_width=2)
                self.all_tris.add(t0, t1)
        
        self.play(FadeIn(self.all_tris), run_time=2)
        self.wait(1)
        self.remove(step4_text)

    def mark_corners(self):
        # Step label
        step5_text = Text("5. Identifying Boundary Vertices", font_size=24, color=RED).next_to(self.title, DOWN)
        self.play(Write(step5_text))
        
        # Corners (0, N-1, N*(N-1), N*N-1) -> (0, 2, 6, 8)
        self.corner_indices = [0, self.N - 1, self.N * (self.N - 1), self.N * self.N - 1]
        
        corner_label = Text("Corners Marked", font_size=20, color=RED).to_edge(RIGHT).shift(UP)
        self.play(Write(corner_label))
        # Animate dots and scale for emphasis
        self.play(*[self.dots[idx].animate.set_color(RED).scale(2.0) for idx in self.corner_indices])
        self.wait(1)
        self.remove(step5_text)

    def mark_edges(self):
        # Inner-loop edges (1, N+0, N+N-1, N+N...) -> (1, 3, 5, 7) for N=3
        edge_indices = []
        for i in range(1, self.N - 1):
            # i, V*i, V*i+N, V*N+i logic maps to inner edges
            edge_indices.extend([i, self.N * i, self.N * i + (self.N - 1), self.N * (self.N - 1) + i])

        edge_label = Text("Edges Iteratively Marked", font_size=20, color=ORANGE).to_edge(RIGHT)
        self.play(Write(edge_label))

        # Iterative marking as requested
        for idx in edge_indices:
            self.play(self.dots[idx].animate.set_color(ORANGE).scale(1.5), run_time=0.5)