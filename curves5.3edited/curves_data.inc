      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
c...Maximum number of atoms in complex
      parameter (n1=30000)

c...Maximum number of subunits
      parameter (n2=2500)

c...Maximum number of nucleotides
      parameter (n3=1000)

c...Delta change of variables for numerical gradient test 
      parameter (deltv=1.d-8)

c...Conversion factor from degrees to radians 
      parameter (cdr=0.017453293d0)	

c...Conversion factor from radians to degrees 
      parameter (crd=57.29577951d0)

c...Constant: Pi
      parameter (pi=3.141592654d0)

c...Constant: Precision for Eigenvalue calculation
      parameter (range=1.d-8)
