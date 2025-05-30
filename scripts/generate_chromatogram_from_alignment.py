# --- Import necessary modules ---
from Bio import SeqIO 
import matplotlib.pyplot as plt
import numpy as np 
from scipy.ndimage import gaussian_filter1d
import matplotlib.patches as patches 
import matplotlib as mpl

# --- Set consistent plotting parameters ---
mpl.rcParams.update({
    'font.family': 'Avenir Next',
    'font.weight': 'demi', 
    'font.size': 8,
    'text.color': '#231F20',
    'axes.labelcolor': '#231F20',
    'xtick.color': '#231F20',
    'ytick.color': '#231F20',
    'axes.edgecolor': '#231F20',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'pdf.fonttype': 42,  
    'text.usetex': False  
})

# --- Define input ABI file (Sanger sequencing trace) ---
filename = 'R29 PPS 12_6358' 

# --- Load ABI trace data ---
record = SeqIO.read(f"./snapgene_alignment_trace/alignment_files/{filename}.ab1", "abi")
trace_data = record.annotations['abif_raw'] 
base_calls = trace_data['PBAS1'].decode()
ploc = trace_data['PLOC1']

# --- Extract electropherogram signal intensities for each nucleotide channel ---
trace_dict = {
    'G': trace_data['DATA9'],  
    'A': trace_data['DATA10'],
    'T': trace_data['DATA11'],
    'C': trace_data['DATA12'],
}

# --- Define channel colors ---
colors = {
    'A': '#80C080',  # green
    'C': '#8080FF',  # blue
    'G': '#404040',  # dark gray
    'T': '#FF8080',  # red
}

# --- Define the target sequence to locate within the chromatogram ---
target_seq = "ccctcgc".upper()      
reference_seq = "ccctcgc".upper()   
target_len = len(target_seq)

# --- Locate the starting index of the target sequence in the called bases ---
target_index = base_calls.find(target_seq)
if target_index == -1:
    raise ValueError("Target sequence not found in base calls.")

# --- Get the chromatogram peak positions for each base in the target sequence ---
peak_positions = ploc[target_index : target_index + target_len]

# --- Define the plotting window around the region of interest (±20 bases) ---
start_crop = peak_positions[0] - 20
end_crop = peak_positions[-1] + 20
x_crop = np.arange(start_crop, end_crop)

fig, ax = plt.subplots(figsize=(1.5, 0.75), dpi=600, facecolor='white')  # Small figure for compact embedding

# --- Plot smoothed chromatogram traces for each nucleotide ---
for base in ['G', 'A', 'C', 'T']:
    raw = np.array(trace_dict[base][start_crop:end_crop])  
    smooth = gaussian_filter1d(raw, sigma=1.0)  
    ax.fill_between(x_crop, 0, smooth, color=colors[base], alpha=0.3)
    ax.plot(x_crop, smooth, color=colors[base], linewidth=1.2)

ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(x_crop[0], x_crop[-1])
ax.set_facecolor('white')
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout(pad=0.5)
plt.savefig(f"snapgene_alignment_trace/output/chromatogram_{filename}.svg", transparent=False)