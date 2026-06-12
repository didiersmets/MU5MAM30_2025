from manim import *

class SquareMeshVisualizer(Scene):
    def construct(self):
        # Configuration matching the C++ code
        N = 3
        V = N + 1
        a_len = 12.0 # Width
        b_len = 8.0 # Height
        
        title = Text("C++ Square Mesh Generation", font_size=40).to_edge(UP)
        self.play(Write(title))

        # Dictionaries to store references to dot and label objects
        dots = {}
        labels = {}

        def get_pos(i, j):
            # Centering the grid using the C++ formula logic
            x = (a_len * i / N) - (a_len / 2)
            y = (b_len * j / N) - (b_len / 2)
            return np.array([x, y, 0])

        ### STEP 1: Vertices and Indices ###
        step1_text = Text("1. Generating Vertices & Indices (Row-Major Map)", font_size=24, color=BLUE).next_to(title, DOWN)
        self.play(Write(step1_text))
        
        # Chronologically iterating exactly as the C++ nested loop does
        for i in range(V):
            for j in range(V):
                idx = V * i + j
                pos = get_pos(i, j)
                
                dot = Dot(pos, radius=0.1, color=WHITE)
                label = Text(str(idx), font_size=20).next_to(dot, DR, buff=0.1)
                
                dots[idx] = dot
                labels[idx] = label
                
                # Adding vertices sequentially
                self.play(FadeIn(dot), Write(label), run_time=0.05)
        
        self.wait(1)

        ### STEP 2: Single Square Triangulation ###
        self.play(Transform(step1_text, Text("2. Triangulating a Single Cell", font_size=24, color=YELLOW).next_to(title, DOWN)))
        
        # Pick a cell in the middle (i=1, j=1)
        focus_i, focus_j = 1, 1
        v = V * focus_i + focus_j
        
        # Highlight the 4 vertices of the current square
        v_idx = [v, v+1, v+1+V, v+V]
        self.play(*[dots[idx].animate.set_color(YELLOW).scale(1.5) for idx in v_idx])
        self.wait(0.5)

        # Triangle 1: v -> v+1 -> v+1+V
        t1_label = Text(f"Tri 1: {v} -> {v+1} -> {v+1+V}", font_size=20, color=BLUE_C).to_edge(LEFT).shift(UP*2)
        self.play(Write(t1_label))
        
        tri1_arrows = VGroup(
            Arrow(dots[v].get_center(), dots[v+1].get_center(), buff=0.1, color=BLUE_C),
            Arrow(dots[v+1].get_center(), dots[v+1+V].get_center(), buff=0.1, color=BLUE_C),
            Arrow(dots[v+1+V].get_center(), dots[v].get_center(), buff=0.1, color=BLUE_C)
        )
        tri1_poly = Polygon(dots[v].get_center(), dots[v+1].get_center(), dots[v+1+V].get_center(), color=BLUE_C, fill_opacity=0.4)
        
        self.play(Create(tri1_arrows), run_time=1.5)
        self.play(FadeIn(tri1_poly))
        self.wait(1)

        # Triangle 2: v -> v+1+V -> v+V
        t2_label = Text(f"Tri 2: {v} -> {v+1+V} -> {v+V}", font_size=20, color=GREEN_C).next_to(t1_label, DOWN, aligned_edge=LEFT)
        self.play(Write(t2_label))

        tri2_arrows = VGroup(
            Arrow(dots[v].get_center(), dots[v+1+V].get_center(), buff=0.1, color=GREEN_C),
            Arrow(dots[v+1+V].get_center(), dots[v+V].get_center(), buff=0.1, color=GREEN_C),
            Arrow(dots[v+V].get_center(), dots[v].get_center(), buff=0.1, color=GREEN_C)
        )
        tri2_poly = Polygon(dots[v].get_center(), dots[v+1+V].get_center(), dots[v+V].get_center(), color=GREEN_C, fill_opacity=0.4)
        
        self.play(Create(tri2_arrows), run_time=1.5)
        self.play(FadeIn(tri2_poly))
        self.wait(1)

        # Clean up arrows and text to show the rest of the mesh
        self.play(FadeOut(tri1_arrows), FadeOut(tri2_arrows), FadeOut(t1_label), FadeOut(t2_label))
        self.play(*[dots[idx].animate.set_color(WHITE).scale(1/1.5) for idx in v_idx])

        ### STEP 3: Complete Triangulation ###
        self.play(Transform(step1_text, Text("3. Triangulating the Remaining Mesh", font_size=24, color=WHITE).next_to(title, DOWN)))
        
        all_tris = VGroup()
        for i in range(N):
            for j in range(N):
                if i == focus_i and j == focus_j:
                    continue # Skip the one we already drew
                
                curr_v = V * i + j
                t1 = Polygon(dots[curr_v].get_center(), dots[curr_v+1].get_center(), dots[curr_v+1+V].get_center(), color=BLUE_E, fill_opacity=0.2, stroke_width=2)
                t2 = Polygon(dots[curr_v].get_center(), dots[curr_v+1+V].get_center(), dots[curr_v+V].get_center(), color=GREEN_E, fill_opacity=0.2, stroke_width=2)
                all_tris.add(t1, t2)
        
        self.play(Create(all_tris), run_time=2)
        self.wait(1)

        ### STEP 4: Boundary Marking ###
        # Dim the triangles to focus on boundaries
        self.play(all_tris.animate.set_opacity(0.05), tri1_poly.animate.set_opacity(0.05), tri2_poly.animate.set_opacity(0.05))
        self.play(Transform(step1_text, Text("4. Marking Boundary Vertices", font_size=24, color=RED).next_to(title, DOWN)))

        # Corners
        corner_indices = [0, N, V*N, V*N + N]
        self.play(*[dots[idx].animate.set_color(RED).scale(1.8) for idx in corner_indices])
        corner_label = Text("Corners Marked", font_size=20, color=RED).to_edge(RIGHT).shift(UP)
        self.play(Write(corner_label))
        self.wait(1)

        # Edges
        edge_indices = []
        for i in range(1, N):
            edge_indices.extend([i, V*i, V*i+N, V*N+i])

        self.play(*[dots[idx].animate.set_color(ORANGE).scale(1.5) for idx in edge_indices])
        edge_label = Text("Edges Marked", font_size=20, color=ORANGE).next_to(corner_label, DOWN)
        self.play(Write(edge_label))
        
        self.wait(3)