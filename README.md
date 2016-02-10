# __float128 sinq(x) precision analyzer
I had some issues with [gcc](https://gcc.gnu.org/) [quadmath's](https://gcc.gnu.org/) sinq(x) on
__float128 types, so I decided to write some code to gather statistics about
the deviation from a reference sin(x) generated with [MPFR](http://www.mpfr.org/). This is the result.

## Usage

Compile and run. Sending HUP will print the stats so far, and CTRL-C will exit and print the final stats.

## Conclusion

The problem must have been in my other project, because the output from this code clearly shows that __float128
and sinq(x) is an improvement in precision over float, double, and long double.
