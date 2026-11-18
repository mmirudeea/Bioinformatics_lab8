import numpy as np
import matplotlib.pyplot as plt

seq_part1 = "CGGACTGATCTATCT"
seq_part2 = "A" * 27
seq_part3 = "CGTAGCATCTATCGATCTATCTAGCATCTA"
S = seq_part1 + seq_part2 + seq_part3

L = 30

def gc_content(window):
    c_count = window.count('C')
    g_count = window.count('G')
    return ((c_count + g_count) / len(window)) * 100

def info_content(window):
    n = len(window)
    if n <= 1: return 0.0
    
    outer_sum = 0.0
    for i in range(n - 1):
        inner_sum = 0.0
        terms_count = n - 1 - i 
        
        for j in range(1, terms_count + 1):
            if window[i] == window[i + j]:
                inner_sum += 1.0
                
        outer_sum += (inner_sum / terms_count) * 100.0
        
    sigma = 100.0 - (outer_sum / (n - 1))
    return sigma

gc_values = []
ic_values = []
num_windows = len(S) - L + 1

for w in range(num_windows):
    window = S[w : w+L]
    gc_values.append(gc_content(window))
    ic_values.append(info_content(window))

GG_mean = np.mean(gc_values)
IG_mean = np.mean(ic_values)
print(f"Overall Average C+G% (GG): {GG_mean:.2f}%")
print(f"Overall Average Info Content (IG): {IG_mean:.2f}")

def center_of_weight(values):
    total_weight = sum(values)
    if total_weight == 0:
        return 0
    weighted_sum = sum(index * val for index, val in enumerate(values))
    return weighted_sum / total_weight

gc_center = center_of_weight(gc_values)
ic_center = center_of_weight(ic_values)

print(f"Center of Weight (C+G%): index {gc_center:.2f}")
print(f"Center of Weight (Info Content): index {ic_center:.2f}")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

x_axis = range(num_windows)
ax1.plot(x_axis, gc_values, label='C+G% Content', color='blue', marker='o', markersize=4)
ax1.plot(x_axis, ic_values, label='Info Content ($\sigma$)', color='red', marker='x', markersize=4)
ax1.set_title("DNA Sliding Window Patterns (L=30)")
ax1.set_xlabel("Window Start Index")
ax1.set_ylabel("Percentage Value")
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.6)

ax2.scatter([gc_center], [1], color='blue', s=100, label=f'C+G% CoW ({gc_center:.2f})')
ax2.scatter([ic_center], [2], color='red', s=100, label=f'Info Content CoW ({ic_center:.2f})')
ax2.axvline(x=gc_center, color='blue', linestyle='--', alpha=0.5)
ax2.axvline(x=ic_center, color='red', linestyle='--', alpha=0.5)

ax2.set_xlim(0, num_windows)
ax2.set_ylim(0, 3)
ax2.set_yticks([1, 2])
ax2.set_yticklabels(['C+G%', 'Info Content'])
ax2.set_title("Center of Weight for Patterns")
ax2.set_xlabel("Window Index Position")
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()