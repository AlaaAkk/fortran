# fortran


To make it work write: bash final.sh then ./main.x.
In the main.f90, I am calling first input_h2.f90 that reads the geometry and basis of H2. 
The size of the variables is given as an input by the user when calling input_h2.f90.
Then overlap_matrix,f90 is called to calculate the overlap matrix. 
After that, the same process for H2O however I call input_h2o.f90 to read the geometry and basis.
