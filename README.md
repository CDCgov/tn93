##We're going to start with a function that we're going to call tn93_distance that takes three arguments - Sequence 1, Sequence 2, and matchMode

##We're also going to need some constants for the calculations. WHich?

 * mapChar - an array with values for each nucleotide (incl ambiguities), the numbers are their order in the resolutions array
 * resolutions - two-dimensional array, the number retrieved from mapChar is the index for the first dimension of this array
 * resolutionsCount - an array with decimal values for the certainty associated with the resolution

##For the distance calculation, the difference is in the matchMode.

 * SKIP is super simple - if they're both unambiguous then add 1 to the appropriate slot in pairwiseCounts
 * GAPMM looks like it turns gaps that are in only one of the two sequencesinto "N"s
 * RESOLVE counts selects the match if possible, averages otherwise
 * AVERAGE just always takes the average of the values

##Critical variables:

 * dist (the tn93 distance)
 * pairwiseCounts (the values used to calculate the distance)

##Algorithm:

Calculate the pairwise nucleotide counts for a sequence pair

 * Get the length of the shortest sequence
 * Step through each location in the sequence
 * Update pairwiseCounts for these two nucleotides
 * Calculate the distance from the pairwiseCounts
  * Get observed nucleotide frequencies
  * Calculate the fraction of the sequence that is not gap
  * Get count of AG mismatches
  * Get count of CT mismatches
  * Calculate tv (don't know what it stands for)
  * If there are any nucleotides that don't appear there's a separate calculation for that
  * Otherwise do the longer calculation
