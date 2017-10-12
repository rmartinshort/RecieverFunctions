These are some test data for the pwaveqn code. The data files
are uln_1995_113_b_057_d_045.[z,n,e]

The output files from pwaveqn are ULN_TEST.[eqr,eqt,aftn]

Here is a script of the deconvolution run. Since I use my own sacio libs,
which don't have setkhv written yet, I get those warning messages.

Script started on Mon Jun 19 08:57:45 2006
[ki:Rftn.Codes/RForward/TestData] cammon% pwaveqn
Specify quake file: uln_1995_113_b_057_d_045.
Real data (y or n)? y
Window data (y or n)? n
 npts= 3001 nft= 4096 fny=10.0000 delf=  0.0049 dt= 0.050
Specify outfil: ULN_TEST
Trough filler, c =  0.0001
Gaussian scale, a = 1.0
Enter phase shift: 30.0
NPTS = 3001
Warning setkhv is not implemented yet.
Warning setkhv is not implemented yet.
Warning setkhv is not implemented yet.
NPTS = 3001
Warning setkhv is not implemented yet.
Warning setkhv is not implemented yet.
Warning setkhv is not implemented yet.
NPTS = 3001
Warning setkhv is not implemented yet.
Warning setkhv is not implemented yet.
Warning setkhv is not implemented yet.
Try another (y or n)? n
[ki:Rftn.Codes/RForward/TestData] cammon% ^D
Script done on Mon Jun 19 08:58:20 2006


Comments on the results: The low-frequency noise in the EQR and EQT suggest
that I should have used a slightly larger water-level. The results suffice
for the test.
