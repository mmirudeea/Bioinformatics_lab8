# Bioinformatics_lab8
## Project Overview
This software application generates DNA patterns from DNA sequence gene parameters. It computes the **C+G% (CpG content)** and the **Information Content** values from each sliding window across the sequence.

## Algorithm & Implementation Steps
To ensure the algorithm implementation works correctly, the following validation steps were followed:

1. **Target Sequence**: Use the DNA sequence `S = CGGACTGATCTATCT(A27)CGTAGCATCTATCGATCTATCTAGCATCTA`.
2. **Sliding Window**: Use a sliding window of length `L = 30b`.
3. **C+G% Calculation**: Build a function to process the CpG content (Target return: `GG = ~29.27`).
4. **Information Content Calculation**: Build a function to process the information content (Target return: `IG = ~27.53`).
5. **Pattern Visualization**: Plot the generated DNA patterns on a chart.
6. **Center of Weight**: Calculate the center of weight for each pattern.
7. **Center of Weight Visualization**: Take the calculated center of each pattern and plot them on a secondary chart.
8. **Verification**: Take the DNA sequence parameter, generate the pattern, and open the `pronk` repository on GitHub to verify if the resulting pattern is identical.
