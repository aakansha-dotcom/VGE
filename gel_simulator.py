import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def calculate_mobility(dna_lengths, control_length=None):
    """
    Calculate mobility with optional control fragment.
    Returns distances in normalized (0-1) and real-world cm (scaled to 8cm gel length).
    """
    all_lengths = np.array(dna_lengths.copy())
    if control_length:
        all_lengths = np.append(all_lengths, control_length)

    mobility_raw = 1.0 / (all_lengths + 1e-6)
    min_mob, max_mob = np.min(mobility_raw), np.max(mobility_raw)

    if max_mob == min_mob:
        norm_dist = np.full_like(mobility_raw, 0.5)
    else:
        norm_dist = (mobility_raw - min_mob) / (max_mob - min_mob)

    # Convert to real-world distances (8cm gel length is standard)
    real_dist_cm = 8 * (1 - norm_dist)  # 0cm at well, 8cm max travel

    if control_length:
        return norm_dist[:-1], real_dist_cm[:-1], norm_dist[-1], real_dist_cm[-1]
    return norm_dist, real_dist_cm, None, None

def visualize_gel(dna_lengths, norm_dist, real_dist_cm,
                 control_norm=None, control_real=None,
                 control_length=None, title="Virtual Gel Electrophoresis"):
    """
    Visualizes gel electrophoresis with:
    - Each fragment aligned with its well
    - Consistent band dimensions
    - Color differentiation
    - Migration distance as position
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Parameters
    num_samples = len(dna_lengths)
    has_control = control_length is not None
    total_lanes = num_samples + (1 if has_control else 0)
    
    # Layout parameters
    lane_width = 0.08
    lane_spacing = 0.9 / total_lanes
    band_height = 0.03
    well_height = 0.02
    well_offset = -0.03
    
    # Color setup - distinct colors for each lane
    colors = plt.cm.tab10(np.linspace(0, 1, total_lanes))
    if has_control:
        colors[0] = [1, 0, 0, 0.8]  # Make control red
    
    # Calculate lane positions
    start_x = 0.05
    lane_positions = [start_x + i*lane_spacing for i in range(total_lanes)]
    
    # Draw wells and bands
    for i in range(total_lanes):
        is_control = has_control and i == 0
        lane_x = lane_positions[i]
        
        # Current fragment data
        if is_control:
            length = control_length
            distance = control_norm
            real_dist = control_real
        else:
            idx = i-1 if has_control else i
            length = dna_lengths[idx]
            distance = norm_dist[idx]
            real_dist = real_dist_cm[idx]
        
        # Draw well (rectangle at top)
        well = Rectangle(
            (lane_x - lane_width/2, well_offset),
            lane_width, well_height,
            facecolor='lightgray', edgecolor='red' if is_control else 'black', lw=1.5
        )
        ax.add_patch(well)
        
        # Draw band (rectangle at migration distance)
        band = Rectangle(
            (lane_x - lane_width/2, distance - band_height/2),
            lane_width, band_height,
            facecolor=colors[i], edgecolor='black', lw=1, alpha=0.8
        )
        ax.add_patch(band)
        
        # Add label
        label_text = f"{length} bp\n{real_dist:.1f} cm"
        ax.text(
            lane_x + lane_width/2 + 0.02, distance,
            label_text,
            verticalalignment='center',
            fontsize=10,
            color=colors[i],
            weight='bold' if is_control else 'normal'
        )
        
        # Add lane number
        ax.text(
            lane_x, well_offset - 0.02,
            f"Lane {i+1}",
            ha='center', va='top', fontsize=9
        )
    
    # Add measurement scale
    scale_x = 0.95
    ax.plot([scale_x, scale_x], [0, 1], color='black', lw=1.5)
    for y in np.linspace(0, 1, 5):
        ax.plot([scale_x, scale_x+0.01], [y, y], color='black', lw=1.5)
        ax.text(scale_x+0.02, y, f"{8*(1-y):.1f} cm", ha='left', va='center', fontsize=10)
    
    # Configure plot
    ax.set_xlim(0, 1)
    ax.set_ylim(well_offset - 0.05, 1.1)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=16, pad=20)
    ax.set_ylabel("Migration Distance", fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Well (0 cm)", "Max (8 cm)"])
    
    # Gel background
    ax.set_facecolor('#FFF8DC')
    plt.grid(False)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("=== Gel Electrophoresis Simulator ===")
    
    # Get sample fragments
    dna_lengths = []
    while True:
        user_input = input("Enter sample DNA length (bp) or 'done': ").strip()
        if user_input.lower() == 'done':
            break
        try:
            length = int(user_input)
            if length <= 0:
                print("Must be positive integer.")
            else:
                dna_lengths.append(length)
        except ValueError:
            print("Invalid input. Enter a number or 'done'.")
    
    if not dna_lengths:
        print("No samples entered. Exiting.")
        exit()
    
    # Get control fragment
    control_length = None
    while True:
        control_input = input("Enter control DNA length (bp) or 'skip': ").strip()
        if control_input.lower() == 'skip':
            break
        try:
            control_length = int(control_input)
            if control_length <= 0:
                print("Must be positive integer.")
            else:
                break
        except ValueError:
            print("Invalid input. Enter a number or 'skip'.")
    
    # Calculate distances
    norm_dist, real_dist_cm, control_norm, control_real = calculate_mobility(
        dna_lengths, control_length)
    
    # Visualize
    title = "Gel Electrophoresis"
    if control_length:
        title += f"\nControl: {control_length} bp (Lane 1)"
    
    visualize_gel(dna_lengths, norm_dist, real_dist_cm,
                 control_norm, control_real, control_length,
                 title=title)