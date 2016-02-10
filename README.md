# __float128 sinq(x) precision analyzer
I had some issues with [gcc](https://gcc.gnu.org/) [quadmath's](https://gcc.gnu.org/) sinq(x) on
__float128 types, so I decided to write some code to gather statistics about
the deviation from a reference sin(x) generated with [MPFR](http://www.mpfr.org/). This is the result.

## Usage

Compile and run. Sending HUP will print the stats so far, and CTRL-C will exit and print the final stats.

## Conclusion

The problem must have been in my other project, because the output from this code clearly shows that __float128
and sinq(x) is an improvement in precision over float, double, and long double.

I had it run for 342677246 rounds, and the mean error and its standard devitation is much smaller for __float128 than for long double. The relative differences are what would be expected for the various bit precisions. Apparently the `+0 <= x < PI/2, uniform` generator craps out eventually. Not sure why.

```
#   Distribution: "+0 <= x < PI/2, non-uniform"   Generator: "float"
Samples: 342677246
Relative difference mean: 1.5613321627e-10
Relative difference variance: 5.9718378264e-17
Relative difference standard deviation: 7.7277667061e-09

#   Distribution: "+0 <= x < PI/2, non-uniform"   Generator: "double"
Samples: 342677246
Relative difference mean: 3.1203133708e-19
Relative difference variance: 4.2857310178e-34
Relative difference standard deviation: 2.0702007192e-17

#   Distribution: "+0 <= x < PI/2, non-uniform"   Generator: "long double"
Samples: 342677246
Relative difference mean: 1.4144321307e-22
Relative difference variance: 1.2624510948e-40
Relative difference standard deviation: 1.1235884900e-20

#   Distribution: "+0 <= x < PI/2, non-uniform"   Generator: "__float128"
Samples: 342677246
Relative difference mean: 2.7505184808e-37
Relative difference variance: 6.8596366826e-70
Relative difference standard deviation: 2.6190908122e-35

#   Distribution: "+0 <= x < PI/2, uniform"   Generator: "float"
Samples: 342677246
Relative difference mean: nan
Relative difference variance: nan
Relative difference standard deviation: nan

#   Distribution: "+0 <= x < PI/2, uniform"   Generator: "double"
Samples: 342677246
Relative difference mean: nan
Relative difference variance: nan
Relative difference standard deviation: nan

#   Distribution: "+0 <= x < PI/2, uniform"   Generator: "long double"
Samples: 342677246
Relative difference mean: nan
Relative difference variance: nan
Relative difference standard deviation: nan

#   Distribution: "+0 <= x < PI/2, uniform"   Generator: "__float128"
Samples: 342677246
Relative difference mean: nan
Relative difference variance: nan
Relative difference standard deviation: nan

#   Distribution: "all floats"   Generator: "float"
Samples: 342677246
Relative difference mean: -8.4791809607e-13
Relative difference variance: 2.9093407176e-16
Relative difference standard deviation: 1.7056789609e-08

#   Distribution: "all floats"   Generator: "double"
Samples: 342677246
Relative difference mean: 5.3650625318e-21
Relative difference variance: 1.1200799189e-33
Relative difference standard deviation: 3.3467595057e-17

#   Distribution: "all floats"   Generator: "long double"
Samples: 342677246
Relative difference mean: -1.1127624593e-24
Relative difference variance: 3.0471923360e-40
Relative difference standard deviation: 1.7456209027e-20

#   Distribution: "all floats"   Generator: "__float128"
Samples: 342677246
Relative difference mean: 6.0357927578e-40
Relative difference variance: 1.1055007082e-69
Relative difference standard deviation: 3.3249070788e-35
```
