NEWS:

09-12-2019: Made the following changes:
1. Added some statistics explanation in the vignette
2. Changed the title of the network visualisation in the draw_network function.

19-11-2019: Made the following changes:
1. Changed the title to be shorter
2. Increased the margin in the network and heatmap figures so they produced nicer looker outputs
3. Changed T/F to TRUE/FALSE (doesn't change functionality)

11-11-2019: Made the following changes:
1. Previously, there would sometimes be the following warning:
"In 1:firstOne : numerical expression has 2 elements: only the first used"
The code still produced the correct output, but it has not been re-written so that 
the warning doesn't appear.

6-11-2019: New submission

This package includes basic functions to build networks of GO terms.  
It has separate functions for building these networks, with and without edge weights.
It also has two functions to help visualize the results with igraph and heatmaps.
