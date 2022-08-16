# Calculating the Tamura-Nei distance for a pair of HIV sequences

Calculating the Tamura-Nei (TN93) distance requires a pair of aligned sequences. Step through the sequence and for each position identify the nucelotide in sequence #1 and the nucleotide in sequence #2 and increment the cell in the pairwise_counts matrix that corresponds to those two nucleotides.

The algorithm also requires a "match mode" that specifies how any ambiguous nucleotides that are encountered are handled.

* SKIP ignores any ambiguous nucelotides
* GAPMM treats gaps appearing in only one sequence as mismatches
* RESOLVE resolves the ambiguity if possible, averages possible values otherwise
* AVERAGE takes the average of possible values

## Algorithm:

Calculate the pairwise nucleotide counts for a sequence pair

 * Get the length of the shortest sequence
 * Step through each location in the sequence
     * Select the current nucleotide from each
     * If the nucleotides are unambiguous, add one to the appropriate pairwise_counts cell
     * If the nucleotides are ambiguous, resolve it according to match_mode
 * Calculate the distance from the pairwise_counts
     * Get observed nucleotide frequencies
     * Calculate the fraction of the sequence that is not gap
     * Get count of AG mismatches
     * Get count of CT mismatches
     * Calculate tv (don't know what it stands for)
     * If there are any nucleotides that don't appear there's a separate calculation for that
     * Otherwise do the longer calculation
