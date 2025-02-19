import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi

input_file = input("Please enter the name of the input CSV file (with extension): ")

data = pd.read_csv(input_file)

print("Available columns:", data.columns)

data.columns = data.columns.str.strip()

if 'Range' not in data.columns or 'Spacer Count' not in data.columns:
    raise KeyError("CSV file must contain 'Range' and 'Spacer Count' columns.")

ranges = data['Range']
spacer_counts = data['Spacer Count']

N = len(ranges)

max_value = max(spacer_counts)
scaling_factor = 0.3
normalized_counts = spacer_counts / max_value * scaling_factor

angles = [(n / float(N) * 2 * pi) + (pi / 2) for n in range(N)]
bar_width = 2 * pi / N * 0.8

fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={'polar': True})

bars = ax.bar(angles, normalized_counts, width=bar_width, bottom=0.5,
              color='brown', edgecolor='black', alpha=0.7)

ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim(0, 1)

ax.spines['polar'].set_visible(False)
ax.grid(False)

output_file = "radial_chart.png"
plt.savefig(output_file, dpi=600, bbox_inches='tight')

plt.tight_layout()
plt.show()
