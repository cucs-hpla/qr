###################################################################
Analysis.md<br />
Eric Peters<br />
cucs-hpla<br />
HW2:Due 2017-03-01<br />
###################################################################

#1.) The starter code does not actually use OpenMP.  Can you annotate loops to improve performance? For which sizes?


#2.) What impediments remain for higher performance?


#3.) A block QR algorithm (known as Tall-Skinny QR, TSQR) is based on the recurrence<br \>

A0 = [Q0 0] [R0] = [Q0 0] [U0] R<br \>
A1 = [0 Q1] [R1] = [0 Q1] [U1] <br \>

I.e., split a tall skinny matrix A into two parts.  Factor each of these parts (can do in parallel).  The result is two right triangular parts, [R0;R1] which can be factored as [U0; U1] R.   This final right triangular matrix is the R corresponding to a Householder QR.  The product of the orthogonal matrices is the Q.  Can you implement a blocked QR algorithm that orthogonalizes a few columns at a time using this method, then applies the block reflector to the remaining panel of the matrix (analagous to Reflect1 in the code)?  What observations can you make about this algorithm?
