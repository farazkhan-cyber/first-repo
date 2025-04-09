import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QLineEdit, QPushButton,
                             QFrame, QGroupBox, QSizePolicy)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QFont, QIcon

# Modern color palette
COLOR_PALETTE = {
    "background": "#2D2D2D",
    "surface": "#404040",
    "primary": "#5680E9",
    "secondary": "#84CEEB",
    "accent": "#8860D0",
    "text": "#E0E0E0",
    "success": "#5ABF95",
    "warning": "#D19A3E",
    "danger": "#BF5A5A"
}

# Configure matplotlib style
plt.style.use('dark_background')
plt.rcParams.update({
    'axes.facecolor': COLOR_PALETTE["surface"],
    'figure.facecolor': COLOR_PALETTE["background"],
    'axes.edgecolor': COLOR_PALETTE["text"],
    'axes.labelcolor': COLOR_PALETTE["text"],
    'text.color': COLOR_PALETTE["text"],
    'xtick.color': COLOR_PALETTE["text"],
    'ytick.color': COLOR_PALETTE["text"],
    'grid.color': '#4A4A4A',
    'lines.linewidth': 2
})

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=12, height=10, dpi=100):
        self.fig = plt.figure(figsize=(width, height), dpi=dpi, facecolor=COLOR_PALETTE["background"])
        gs = GridSpec(5, 2, figure=self.fig,
                      height_ratios=[1, 1, 1, 1.5, 1.5],
                      width_ratios=[3, 1])
        
        # Create the styled axes
        self.axs = {
            'beam': self._styled_axis(gs[0, 0]),
            'reaction_A': self._styled_axis(gs[1, 0]),
            'reaction_B': self._styled_axis(gs[2, 0]),
            'shear': self._styled_axis(gs[3, 0]),
            'moment': self._styled_axis(gs[4, 0]),
            'legend': self._styled_axis(gs[0:2, 1], grid=False),
            'envelope_sf': self._styled_axis(gs[3, 1]),
            'envelope_bm': self._styled_axis(gs[4, 1])
        }
        
        super().__init__(self.fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.updateGeometry()
        plt.tight_layout()

    def _styled_axis(self, position, grid=True):
        ax = self.fig.add_subplot(position)
        ax.set_facecolor(COLOR_PALETTE["surface"])
        for spine in ax.spines.values():
            spine.set_color(COLOR_PALETTE["text"])
        ax.tick_params(colors=COLOR_PALETTE["text"])
        ax.xaxis.label.set_color(COLOR_PALETTE["text"])
        ax.yaxis.label.set_color(COLOR_PALETTE["text"])
        ax.title.set_color(COLOR_PALETTE["primary"])
        if grid:
            ax.grid(True, color='#4A4A4A', linestyle='--')
        return ax

class BeamAnalysisApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Advanced Beam Analysis System")
        self.setGeometry(100, 100, 1400, 1000)
        # Load window icon with proper path handling

        path = os.path.dirname(os.path.realpath(__file__))+"/resources/logo.png"
        if os.path.exists(path):
            self.setWindowIcon(QIcon(path))
        else:
            print(f"Warning: Icon not found") 
        # Animation control
        self.anim = None
        self.is_animating = False
        self.max_values = {
            'SF': {'value': 0, 'position': 0},
            'BM': {'value': 0, 'position': 0}
        }
        
        # Create the main central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout (horizontal)
        main_layout = QHBoxLayout(central_widget)
        
        # Create GUI components
        self.create_input_panel(main_layout)
        self.create_plots(main_layout)
        
    def create_input_panel(self, main_layout):
        # Create left panel for inputs and results
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(20, 20, 20, 20)
        
        # Input parameters
        input_group = QGroupBox("Input Parameters")
        input_layout = QVBoxLayout(input_group)
        input_layout.setContentsMargins(15, 15, 15, 15)
        
        self.entries = {}
        fields = [
            ('Beam Length (L, m)', 'L', '10'),
            ('Load W1 (kN)', 'W1', '50'),
            ('Load W2 (kN)', 'W2', '30'),
            ('Load Spacing (x, m)', 'x', '3'),
            ('Speed (ms/frame)', 'speed', '40')
        ]
        
        for label, key, default in fields:
            row_layout = QHBoxLayout()
            label_widget = QLabel(label)
            label_widget.setFixedWidth(160)
            row_layout.addWidget(label_widget)
            
            self.entries[key] = QLineEdit()
            self.entries[key].setText(default)
            row_layout.addWidget(self.entries[key])
            
            input_layout.addLayout(row_layout)
        
        left_layout.addWidget(input_group)
        
        # Control buttons
        btn_frame = QFrame()
        btn_layout = QHBoxLayout(btn_frame)
        btn_layout.setContentsMargins(0, 15, 0, 15)
        
        start_button = QPushButton("Start")
        start_button.clicked.connect(self.start_animation)
        btn_layout.addWidget(start_button)
        
        stop_button = QPushButton("Stop")
        stop_button.clicked.connect(self.stop_animation)
        btn_layout.addWidget(stop_button)
        
        reset_button = QPushButton("Reset")
        reset_button.clicked.connect(self.reset_analysis)
        btn_layout.addWidget(reset_button)
        
        left_layout.addWidget(btn_frame)
        
        # Results display
        result_group = QGroupBox("Maximum Values")
        result_layout = QVBoxLayout(result_group)
        result_layout.setContentsMargins(15, 15, 15, 15)
        
        self.result_labels = {
            'max_SF': QLabel("Max Shear: -"),
            'max_BM': QLabel("Max Moment: -"),
            'current_SF': QLabel("Current Shear: -"),
            'current_BM': QLabel("Current Moment: -"),
            'current_position': QLabel("Current Position: -")
        }
        
        for lbl in self.result_labels.values():
            result_layout.addWidget(lbl)
            
        left_layout.addWidget(result_group)
        
        # Information panel
        info_group = QGroupBox("Analysis Information")
        info_layout = QVBoxLayout(info_group)
        info_layout.setContentsMargins(15, 15, 15, 15)
        
        info_content = """This application simulates beam behavior under moving loads.

Color Legend:
• Red square: First load (W₁)
• Blue square: Second load (W₂)
• Green line: Shear force diagram
• Purple line: Bending moment diagram

The simulation tracks maximum values and their positions along the beam during the entire load movement cycle."""
        
        info_label = QLabel(info_content)
        info_label.setWordWrap(True)
        info_layout.addWidget(info_label)
        
        left_layout.addWidget(info_group)
        
        # Add left panel to main layout with proper sizing
        left_panel.setMaximumWidth(340)
        main_layout.addWidget(left_panel)
        
    def create_plots(self, main_layout):
        # Right panel for plots
        self.canvas = MplCanvas(self)
        self.initialize_plots()
        
        # Add the canvas to the main layout
        main_layout.addWidget(self.canvas, 1)
        
    def initialize_plots(self):
        # Get the axes objects from the canvas
        self.axs = self.canvas.axs
        
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
        
        # Initialize annotations with contrasting text color for visibility
        self.shear_max_text = self.axs['shear'].text(
            0.05, 0.90, '', transform=self.axs['shear'].transAxes,
            fontsize=9, color='black', bbox=dict(facecolor='white', alpha=0.8))
        self.moment_max_text = self.axs['moment'].text(
            0.05, 0.90, '', transform=self.axs['moment'].transAxes,
            fontsize=9, color='black', bbox=dict(facecolor='white', alpha=0.8))
        
        # Setup envelope diagrams
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
        
        # Legend/Info panel
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
        
        self.canvas.fig.tight_layout()
        
    def validate_inputs(self):
        try:
            return {
                'L': max(2, float(self.entries['L'].text())),
                'W1': max(1, float(self.entries['W1'].text())),
                'W2': max(1, float(self.entries['W2'].text())),
                'x': max(1, float(self.entries['x'].text())),
                'speed': max(20, int(self.entries['speed'].text()))
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
        for key, lbl in self.result_labels.items():
            label_text = lbl.text().split(':')[0]
            lbl.setText(f"{label_text}: -")
            
    def start_animation(self):
        if self.is_animating: 
            return
        
        params = self.validate_inputs()
        if not params: 
            return
        
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
            self.canvas.fig, update,
            frames=np.linspace(0, params['L'] - params['x'], 100),
            interval=params['speed'],
            blit=True,
            repeat=True,
            cache_frame_data=False
        )
        
        self.canvas.draw()
        
    def update_result_labels(self, current_pos, current_sf, current_bm):
        # Update maximum value displays
        self.result_labels['max_SF'].setText(
            f"Max Shear: {self.max_values['SF']['value']:.1f} kN @ {self.max_values['SF']['position']:.1f}m")
        self.result_labels['max_BM'].setText(
            f"Max Moment: {self.max_values['BM']['value']:.1f} kN·m @ {self.max_values['BM']['position']:.1f}m")
        
        # Update current value displays
        self.result_labels['current_SF'].setText(f"Current Shear: {current_sf:.1f} kN")
        self.result_labels['current_BM'].setText(f"Current Moment: {current_bm:.1f} kN·m")
        self.result_labels['current_position'].setText(f"Current Position: {current_pos:.1f} m")
        
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
    app = QApplication(sys.argv)
    window = BeamAnalysisApp()
    window.show()
    sys.exit(app.exec())
