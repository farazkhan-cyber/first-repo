import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches

class BeamAnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Advanced Beam Analysis System")
        self.root.geometry("1400x1000")
        
        # Animation control
        self.anim = None
        self.is_animating = False
        self.max_values = {
            'SF': {'value': 0, 'position': 0},
            'BM': {'value': 0, 'position': 0}
        }
        
        # Configure root grid to expand properly
        self.root.columnconfigure(0, weight=1)  # Input panel
        self.root.columnconfigure(1, weight=4)  # Plot panel
        self.root.rowconfigure(0, weight=1)
        
        # Create GUI components
        self.create_input_panel()
        self.create_plots()
        
    def create_input_panel(self):
        # Left panel for inputs and results
        left_frame = ttk.Frame(self.root, padding=20)
        left_frame.grid(row=0, column=0, sticky='nsew')
        
        # Input parameters
        input_frame = ttk.LabelFrame(left_frame, text="Input Parameters", padding=15)
        input_frame.pack(fill='x', pady=10)
        
        self.entries = {}
        fields = [
            ('Beam Length (L, m)', 'L', '10'),
            ('Load W1 (kN)', 'W1', '50'),
            ('Load W2 (kN)', 'W2', '30'),
            ('Load Spacing (x, m)', 'x', '3'),
            ('Speed (ms/frame)', 'speed', '40')
        ]
        
        for label, key, default in fields:
            row = ttk.Frame(input_frame)
            row.pack(fill='x', pady=4)
            ttk.Label(row, text=label, width=22).pack(side='left')
            self.entries[key] = ttk.Entry(row)
            self.entries[key].pack(side='right', expand=True, fill='x')
            self.entries[key].insert(0, default)
        
        # Control buttons
        btn_frame = ttk.Frame(left_frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Start", command=self.start_animation).pack(side='left', padx=7)
        ttk.Button(btn_frame, text="Stop", command=self.stop_animation).pack(side='left', padx=7)
        ttk.Button(btn_frame, text="Reset", command=self.reset_analysis).pack(side='left', padx=7)
        
        # Results display
        result_frame = ttk.LabelFrame(left_frame, text="Maximum Values", padding=15)
        result_frame.pack(fill='both', expand=True, pady=10)
        self.result_labels = {
            'max_SF': ttk.Label(result_frame, text="Max Shear: -"),
            'max_BM': ttk.Label(result_frame, text="Max Moment: -"),
            'current_SF': ttk.Label(result_frame, text="Current Shear: -"),
            'current_BM': ttk.Label(result_frame, text="Current Moment: -"),
            'current_position': ttk.Label(result_frame, text="Current Position: -")
        }
        for lbl in self.result_labels.values():
            lbl.pack(anchor='w', pady=4)
            
        # Information panel - IMPROVED LAYOUT
        info_frame = ttk.LabelFrame(left_frame, text="Analysis Information", padding=15)
        info_frame.pack(fill='both', expand=True, pady=10)
        
        info_content = """This application simulates beam behavior under moving loads.

Color Legend:
• Red square: First load (W1)
• Blue square: Second load (W2)
• Green line: Shear force diagram
• Purple line: Bending moment diagram

The simulation tracks maximum values and their positions along the beam during the entire load movement cycle."""
        
        ttk.Label(info_frame, text=info_content, wraplength=300, justify='left').pack(fill='both', expand=True)
        
    def create_plots(self):
        # Right panel for plots
        plot_frame = ttk.Frame(self.root)
        plot_frame.grid(row=0, column=1, sticky='nsew')
        
        # Configure figure with enhanced layout - IMPROVED GRIDSPEC LAYOUT
        self.fig = plt.figure(figsize=(12, 10), dpi=100)
        gs = GridSpec(5, 2, figure=self.fig, 
                     height_ratios=[1, 1, 1, 1.5, 1.5],
                     width_ratios=[3, 1])  # Main plots and side panel
        
        # Create subplots with adjusted positioning
        self.axs = {
            'beam': self.fig.add_subplot(gs[0, 0]),
            'reaction_A': self.fig.add_subplot(gs[1, 0]),
            'reaction_B': self.fig.add_subplot(gs[2, 0]),
            'shear': self.fig.add_subplot(gs[3, 0]),
            'moment': self.fig.add_subplot(gs[4, 0]),
            'legend': self.fig.add_subplot(gs[0:2, 1]),  # EXPANDED LEGEND SPACE
            'envelope_sf': self.fig.add_subplot(gs[3, 1]),  # ALIGNED WITH SHEAR GRAPH
            'envelope_bm': self.fig.add_subplot(gs[4, 1])   # ALIGNED WITH MOMENT GRAPH
        }
        
        # Initialize plots
        self.initialize_plots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        
    def initialize_plots(self):
        # Beam diagram
        self.axs['beam'].set_title("Moving Loads on Beam", fontsize=11, pad=10)
        self.beam_line, = self.axs['beam'].plot([], [], 'k-', lw=3)
        
        # Support markers
        self.support_A = self.axs['beam'].scatter([], [], marker='^', s=100, color='black')
        self.support_B = self.axs['beam'].scatter([], [], marker='^', s=100, color='black')
        
        # Loads
        self.load1, = self.axs['beam'].plot([], [], 's', markersize=12, markerfacecolor='red', markeredgecolor='black')
        self.load2, = self.axs['beam'].plot([], [], 's', markersize=12, markerfacecolor='blue', markeredgecolor='black')
        
        # Annotations
        self.beam_length_text = self.axs['beam'].text(0, 0, '', ha='center')
        self.load_spacing_text = self.axs['beam'].text(0, 0, '', ha='center')
        self.load1_text = self.axs['beam'].text(0, 0, '', ha='center')
        self.load2_text = self.axs['beam'].text(0, 0, '', ha='center')
        self.reaction_A_text = self.axs['beam'].text(0, 0, '', ha='center')
        self.reaction_B_text = self.axs['beam'].text(0, 0, '', ha='center')
        
        # Turn off axis
        self.axs['beam'].axis('off')
        
        # Reaction plots
        for reaction in ['reaction_A', 'reaction_B']:
            self.axs[reaction].clear()
            self.axs[reaction].grid(True, linestyle=':', alpha=0.6)
            self.axs[reaction].set_ylabel('Reaction (kN)', fontsize=9)
            
        self.axs['reaction_A'].set_title("Reaction at A", fontsize=10)
        self.axs['reaction_B'].set_title("Reaction at B", fontsize=10)
        self.ra_line, = self.axs['reaction_A'].plot([], [], 'r-', lw=1.5)
        self.rb_line, = self.axs['reaction_B'].plot([], [], 'b-', lw=1.5)
        
        # Current position marker for reactions
        self.ra_marker, = self.axs['reaction_A'].plot([], [], 'ro', markersize=8)
        self.rb_marker, = self.axs['reaction_B'].plot([], [], 'bo', markersize=8)
        
        # Shear and moment diagrams
        for diagram in ['shear', 'moment']:
            self.axs[diagram].clear()
            self.axs[diagram].grid(True, linestyle=':', alpha=0.6)
            self.axs[diagram].set_xlabel('Beam Position (m)', fontsize=9)
            
        self.axs['shear'].set_title("Shear Force Diagram", fontsize=10)
        self.axs['moment'].set_title("Bending Moment Diagram", fontsize=10)
        self.shear_line, = self.axs['shear'].plot([], [], 'g-', lw=2)
        self.moment_line, = self.axs['moment'].plot([], [], 'm-', lw=2)
        
        # Initialize annotations
        self.shear_max_text = self.axs['shear'].text(0.05, 0.90, '', transform=self.axs['shear'].transAxes,
                                                    fontsize=9, bbox=dict(facecolor='white', alpha=0.8))
        self.moment_max_text = self.axs['moment'].text(0.05, 0.90, '', transform=self.axs['moment'].transAxes,
                                                      fontsize=9, bbox=dict(facecolor='white', alpha=0.8))
        
        # Setup envelope diagrams - IMPROVED TITLES
        self.axs['envelope_sf'].clear()
        self.axs['envelope_sf'].set_title("SF Envelope", fontsize=10)
        self.axs['envelope_sf'].grid(True, linestyle=':', alpha=0.6)
        self.max_sf_envelope, = self.axs['envelope_sf'].plot([], [], 'g-', lw=2, label="Max")
        self.min_sf_envelope, = self.axs['envelope_sf'].plot([], [], 'g--', lw=2, label="Min")
        self.axs['envelope_sf'].legend(loc='upper right', fontsize=8)
        
        self.axs['envelope_bm'].clear()
        self.axs['envelope_bm'].set_title("BM Envelope", fontsize=10)
        self.axs['envelope_bm'].grid(True, linestyle=':', alpha=0.6)
        self.axs['envelope_bm'].set_xlabel('Beam Position (m)', fontsize=9)
        self.max_bm_envelope, = self.axs['envelope_bm'].plot([], [], 'm-', lw=2, label="Max")
        self.axs['envelope_bm'].legend(loc='upper right', fontsize=8)
        
        # IMPROVED LEGEND/INFO PANEL with better spacing and organization
        self.axs['legend'].clear()
        self.axs['legend'].set_title("Beam Analysis Guide", fontsize=11)
        self.axs['legend'].axis('off')
        
        # Add visual legend with improved spacing
        self.axs['legend'].text(0.05, 0.95, "Legend:", fontsize=10, fontweight='bold')
        self.axs['legend'].text(0.1, 0.88, "● Red Load (W₁)", color='red', fontsize=9)
        self.axs['legend'].text(0.1, 0.81, "● Blue Load (W₂)", color='blue', fontsize=9)
        self.axs['legend'].text(0.1, 0.74, "● Green: Shear", color='green', fontsize=9)
        self.axs['legend'].text(0.1, 0.67, "● Purple: Moment", color='purple', fontsize=9)
        
        # Add formulas with improved spacing and clarity
        self.axs['legend'].text(0.05, 0.55, "Key Formulas:", fontsize=10, fontweight='bold')
        
        # Reaction formulas with better spacing
        self.axs['legend'].text(0.1, 0.48, "Reactions:", fontsize=9, fontstyle='italic')
        self.axs['legend'].text(0.15, 0.41, "Ra = W₁(1-a/L) + W₂(1-b/L)", fontsize=8)
        self.axs['legend'].text(0.15, 0.34, "Rb = W₁(a/L) + W₂(b/L)", fontsize=8)
        
        # Moment formulas with better spacing
        self.axs['legend'].text(0.1, 0.27, "Moment:", fontsize=9, fontstyle='italic')
        self.axs['legend'].text(0.15, 0.20, "M(x) = Ra·x - W₁·max(0,x-a)", fontsize=8)
        self.axs['legend'].text(0.15, 0.13, "        - W₂·max(0,x-b)", fontsize=8)
        
        # Add variables explanation
        self.axs['legend'].text(0.1, 0.06, "Where a, b = load positions", fontsize=8)
        
        # Setup envelopes storage
        self.envelope_data = {
            'x': [],
            'max_sf': [],
            'min_sf': [],
            'max_bm': []
        }
        
        plt.tight_layout()
        
    def validate_inputs(self):
        try:
            return {
                'L': max(2, float(self.entries['L'].get())),
                'W1': max(1, float(self.entries['W1'].get())),
                'W2': max(1, float(self.entries['W2'].get())),
                'x': max(1, float(self.entries['x'].get())),
                'speed': max(20, int(self.entries['speed'].get()))
            }
        except ValueError:
            return None
        
    def reset_analysis(self):
        self.stop_animation()
        for ax in self.axs.values():
            ax.clear()
        self.initialize_plots()
        self.canvas.draw()
        
        # Reset result labels
        for lbl in self.result_labels.values():
            lbl.config(text=lbl['text'].split(':')[0] + ': -')
            
    def start_animation(self):
        if self.is_animating: return
        
        params = self.validate_inputs()
        if not params: return
        
        # Ensure load spacing is valid
        params['x'] = min(params['x'], params['L'] - 1)
        
        # Reset system
        self.max_values = {'SF': {'value': 0, 'position': 0},
                          'BM': {'value': 0, 'position': 0}}
        
        # Reset envelope data
        self.envelope_data = {
            'x': np.linspace(0, params['L'], 200),
            'max_sf': np.zeros(200),
            'min_sf': np.zeros(200),
            'max_bm': np.zeros(200)
        }
        
        # Animation variables
        self.is_animating = True
        self.history = {'pos': [], 'Ra': [], 'Rb': []}
        
        def update(frame):
            current_pos = frame % (params['L'] - params['x'])
            a = current_pos
            b = current_pos + params['x']
            
            # Ensure loads stay on beam
            if b > params['L']:
                a = params['L'] - params['x']
                b = params['L']
            
            # Calculate reactions
            Ra = params['W1']*(1 - a/params['L']) + params['W2']*(1 - b/params['L'])
            Rb = params['W1']*(a/params['L']) + params['W2']*(b/params['L'])
            
            # Update beam visualization
            self.axs['beam'].set_xlim(-1, params['L'] + 1)
            self.axs['beam'].set_ylim(-3, 3)
            
            # Draw beam and supports
            self.beam_line.set_data([0, params['L']], [0, 0])
            self.support_A.set_offsets([[0, -0.2]])
            self.support_B.set_offsets([[params['L'], -0.2]])
            
            # Draw loads
            self.load1.set_data([a], [0])
            self.load2.set_data([b], [0])
            
            # Add annotations
            # Beam length annotation
            self.beam_length_text.set_position((params['L']/2, -1.2))
            self.beam_length_text.set_text(f"L = {params['L']} m")
            
            # Load spacing annotation
            self.load_spacing_text.set_position(((a + b)/2, 0.5))
            self.load_spacing_text.set_text(f"{params['x']} m")
            
            # Load values
            self.load1_text.set_position((a, 1.5))
            self.load1_text.set_text(f"{params['W1']} kN")
            
            self.load2_text.set_position((b, 1.5))
            self.load2_text.set_text(f"{params['W2']} kN")
            
            # Reaction values
            self.reaction_A_text.set_position((0, -2))
            self.reaction_A_text.set_text(f"Ra = {Ra:.1f} kN")
            
            self.reaction_B_text.set_position((params['L'], -2))
            self.reaction_B_text.set_text(f"Rb = {Rb:.1f} kN")
            
            # Track maximum values
            x_vals = np.linspace(0, params['L'], 200)  # More points for smoother curves
            shear, moment = [], []
            current_max_sf = 0
            current_max_bm = 0
            max_sf_pos = 0
            max_bm_pos = 0
            
            for i, xi in enumerate(x_vals):
                # Shear calculation
                if xi < a:
                    V = Ra
                elif a <= xi < b:
                    V = Ra - params['W1']
                else:
                    V = Ra - params['W1'] - params['W2']
                shear.append(V)
                
                # Update envelope for shear force
                self.envelope_data['max_sf'][i] = max(self.envelope_data['max_sf'][i], V)
                self.envelope_data['min_sf'][i] = min(self.envelope_data['min_sf'][i], V)
                
                # Correct moment calculation
                if xi < a:
                    M = Ra * xi
                elif a <= xi < b:
                    M = Ra * xi - params['W1'] * (xi - a)
                else:
                    M = Ra * xi - params['W1'] * (xi - a) - params['W2'] * (xi - b)
                moment.append(M)
                
                # Update envelope for bending moment
                self.envelope_data['max_bm'][i] = max(self.envelope_data['max_bm'][i], M)
                
                # Track maxima
                if abs(V) > current_max_sf:
                    current_max_sf = abs(V)
                    max_sf_pos = xi
                if M > current_max_bm:
                    current_max_bm = M
                    max_bm_pos = xi
            
            # Update global maxima
            if current_max_sf > self.max_values['SF']['value']:
                self.max_values['SF'].update({'value': current_max_sf, 'position': max_sf_pos})
            if current_max_bm > self.max_values['BM']['value']:
                self.max_values['BM'].update({'value': current_max_bm, 'position': max_bm_pos})
            
            # Update results display
            self.update_result_labels(current_pos, current_max_sf, current_max_bm)
            
            # Update reaction plots
            self.update_reaction_plots(current_pos, Ra, Rb)
            
            # Update diagrams
            self.update_diagrams(x_vals, shear, moment)
            
            # Update envelopes
            self.update_envelopes()
            
            return (self.beam_line, self.load1, self.load2, self.ra_line, self.rb_line,
                    self.shear_line, self.moment_line, self.ra_marker, self.rb_marker,
                    self.beam_length_text, self.load_spacing_text, self.load1_text, 
                    self.load2_text, self.reaction_A_text, self.reaction_B_text,
                    self.max_sf_envelope, self.min_sf_envelope, self.max_bm_envelope)
        
        # Create animation
        self.anim = FuncAnimation(
            self.fig, update,
            frames=np.linspace(0, params['L'] - params['x'], 100),
            interval=params['speed'],
            blit=True,
            repeat=True,
            cache_frame_data=False
        )
        
        self.canvas.draw()
        
    def update_result_labels(self, current_pos, current_sf, current_bm):
        # Update maximum value displays
        self.result_labels['max_SF'].config(
            text=f"Max Shear: {self.max_values['SF']['value']:.1f} kN @ {self.max_values['SF']['position']:.1f}m")
        self.result_labels['max_BM'].config(
            text=f"Max Moment: {self.max_values['BM']['value']:.1f} kN·m @ {self.max_values['BM']['position']:.1f}m")
        
        # Update current value displays
        self.result_labels['current_SF'].config(text=f"Current Shear: {current_sf:.1f} kN")
        self.result_labels['current_BM'].config(text=f"Current Moment: {current_bm:.1f} kN·m")
        self.result_labels['current_position'].config(text=f"Current Position: {current_pos:.1f} m")
        
        # Update plot annotations
        self.shear_max_text.set_text(
            f"Max: {self.max_values['SF']['value']:.1f} kN\n@ {self.max_values['SF']['position']:.1f}m")
        self.moment_max_text.set_text(
            f"Max: {self.max_values['BM']['value']:.1f} kN·m\n@ {self.max_values['BM']['position']:.1f}m")
        
    def update_reaction_plots(self, pos, Ra, Rb):
        # Update reaction histories
        self.history['pos'].append(pos)
        self.history['Ra'].append(Ra)
        self.history['Rb'].append(Rb)
        
        # Keep only last 100 points
        if len(self.history['pos']) > 100:
            for key in self.history:
                self.history[key].pop(0)
        
        # Update reaction plots
        self.ra_line.set_data(self.history['pos'], self.history['Ra'])
        self.rb_line.set_data(self.history['pos'], self.history['Rb'])
        
        # Update current position markers
        self.ra_marker.set_data([pos], [Ra])
        self.rb_marker.set_data([pos], [Rb])
        
        # Auto-scale reaction plots
        for reaction in ['reaction_A', 'reaction_B']:
            self.axs[reaction].relim()
            self.axs[reaction].autoscale_view()
            
            # Set x-axis limits to match animation range
            params = self.validate_inputs()
            if params:
                self.axs[reaction].set_xlim(0, params['L'] - params['x'])
            
    def update_diagrams(self, x_vals, shear, moment):
        # Update diagram data
        self.shear_line.set_data(x_vals, shear)
        self.moment_line.set_data(x_vals, moment)
        
        # Set dynamic limits with 20% buffer
        sf_buffer = max(1, self.max_values['SF']['value'] * 0.2)
        bm_buffer = max(1, self.max_values['BM']['value'] * 0.2)
        
        self.axs['shear'].set_ylim(-(self.max_values['SF']['value'] + sf_buffer), 
                                  self.max_values['SF']['value'] + sf_buffer)
        self.axs['moment'].set_ylim(-bm_buffer, 
                                   self.max_values['BM']['value'] + bm_buffer)
        
        # Set x-axis limits to match beam length
        for diagram in ['shear', 'moment']:
            self.axs[diagram].set_xlim(0, x_vals[-1])
            
    def update_envelopes(self):
        # Update envelope plots
        self.max_sf_envelope.set_data(self.envelope_data['x'], self.envelope_data['max_sf'])
        self.min_sf_envelope.set_data(self.envelope_data['x'], self.envelope_data['min_sf'])
        self.max_bm_envelope.set_data(self.envelope_data['x'], self.envelope_data['max_bm'])
        
        # Set y-axis limits for envelope plots
        params = self.validate_inputs()
        if params:
            # Shear envelope limits
            max_sf = np.max(self.envelope_data['max_sf'])
            min_sf = np.min(self.envelope_data['min_sf'])
            sf_range = max(abs(max_sf), abs(min_sf))
            self.axs['envelope_sf'].set_ylim(-sf_range*1.2, sf_range*1.2)
            self.axs['envelope_sf'].set_xlim(0, params['L'])
            
            # Moment envelope limits
            max_bm = np.max(self.envelope_data['max_bm'])
            self.axs['envelope_bm'].set_ylim(-max_bm*0.1, max_bm*1.2)
            self.axs['envelope_bm'].set_xlim(0, params['L'])
        
    def stop_animation(self):
        if self.anim:
            self.anim.event_source.stop()
        self.is_animating = False

if __name__ == "__main__":
    root = tk.Tk()
    app = BeamAnalysisApp(root)
    root.mainloop()

