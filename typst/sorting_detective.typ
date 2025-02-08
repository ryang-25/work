= Sorting Detective

Roland Yang


#table(
  columns: (1fr, 1fr, 1fr),
  inset: 6pt,
  align: (horizon, horizon, center),
  table.header(
    [*Sort*], [*Name*], [*Observed $O(n)$*],
  ),
  [ Alpha ],
  [ Selection ],
  $ O(n^2) $,

  [ Beta ],
  [ Insertion ],
  $ O(n^2) $,

  [ Gamma ],
  [ Bubble ],
  $ O(n^2) $,

  [ Delta ],
  [ Shell ],
  [ $O(n log n)$ or $O(n log^2 n)$ ],

  [ Epsilon ],
  [ Merge ],
  [ $ O(n log n) $ ],
  
  [ Zeta ],
  [ Quick ],
  [ $ O (n log n) $ ],
  
  [ Theta ],
  [ Heap ],
  [ $ O (n log n) $ ],
)


== Methodology

I ended up importing `java.util.Random` and creating a `int[][]` of random
numbers, where the arrays ranged in size from 100-12800 elements. I tested for
20000 runs and measured the time with nanosecond precision using
`System.nanoTime()`. The first 10000 runs were discarded as they were used to
warm up the just-in-time (JIT) compiler to the fastest native tier 4 code. The
code can be found at the end of the report.

#image("assets/sheet1.png")
#align(center)[_SortDetective.java Results_]

For the detective portion of the project, I ended up simply counting the
comparisons to determine the $O(n)$ time complexity. For the most part, the
sorts of $O(n log n)$ time performed similarly, almost indistinguishable from
each other on the graph. Only when extrapolating the trend line to $>=100000$
elements is there an observable difference. Quicksort was the fastest, and the
$O(n^2)$ time complexity sorts were similar to the first graph.

#image("assets/sheet2.png")


== Addendum: SortTester.java

```java
/*
 * Project title: Etude01: Sorting Detective: Chasing Algorithms
 * Roland Yang
 * On my honor, I have neither given nor received unauthorized help on this
 * assignment.
 */

import java.util.Random;

public class SortTester {
    public static void main(String[] args) {

        // stream of randomness
        var rand = new Random();
        var problem = new SortProblem();

        int[][] nums = {
            rand.ints(100).toArray(), rand.ints(200).toArray(),
            rand.ints(400).toArray(), rand.ints(800).toArray(),
            rand.ints(1600).toArray(), rand.ints(3200).toArray(),
            rand.ints(6400).toArray(), rand.ints(12800).toArray()
        };

        int runs = 10;
        // Run for enough iterations to hit tier 4 JIT compile
        for (int i = 0; i < nums.length; i++) {
            long start = 0;
            for (int j = -10000; j < runs; j++) {
                // new array for each iteration
                var a = nums[i].clone();
                // start nanosecond timing
                if (j == 0) start = System.nanoTime();
                switch (args[0]) {
                case "0":
                    problem.selectionSort(a);
                    break;
                case "1":
                    problem.bubbleSort(a);
                    break;
                case "2":
                    problem.insertionSort(a);
                    break;
                case "3":
                    problem.mergeSort(a);
                    break;
                case "4":
                default:
                    problem.quickSort(a);
                    break;
                }
            }
            System.out.print((System.nanoTime() - start)/runs);
            System.out.print(" ");
        }

    }
}
```
