
from bls import load_magnetic_modes, first_order_intensities

# Load the magnon dynamics calculated in Mathematica
data = load_magnetic_modes("results/modes/bloch_gamma.json")

# Calculate first-order BLS intensities
# Incoming polarization: y
# Outgoing polarization: x
weights = first_order_intensities(
    data,
    polarization_in=[0, 1, 0],
    polarization_out=[1, 0, 0],
)

# Print mode frequencies and relative intensities
for i, (frequency, intensity) in enumerate(
    zip(data.frequencies, weights)
):
    print(
        f"Mode {i}: "
        f"frequency = {frequency:.4f}, "
        f"intensity = {intensity:.6f}"
    )

# Check the exported normalization
print("\nBdG norms:", data.bdg_norms)
print("Legacy norms:", data.legacy_norms)