These are individual tests intended to check for reasonable numeric results, independent of comparisons with Acoustic Toolbox BELLHOP.

eigen_arrival: Test the 'V' and 'v' modes that generate both arrivals and associated eigenrays. Mostly a test that the output files are generated and parse. compare the files from 'V', 'E', and 'A' modes (and/or associated binary versions with lower case letters). Note that it works around bugs in how influence is identified by using two receivers.

analysis_3d: Trivial case test that arrivals for receivers along a straight line with constant sound speed follow distance equals rate times time. The third column of the output .arr file should increase linearly. See example output in the folder (example_analysis_3d.arr.txt)
