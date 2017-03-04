## OpenMP Annotations

The only loop I could parallelize for an performance benefit was applying the elementary reflector. The
following chart shows M/N vs speedup from 4 threads for a variety of M/N ratios. I made sure M and N
were large enough that it took at least a second to solve the problem. This should ensure there is
enough work.

![](size_ratio.png)
