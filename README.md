# f2matlab

README file for f2matlab.m

CONTENTS:
-1. SUPPORT f2matlab AND CONSULTING
 0. DISCLAIMER
 1. OBJECTIVE
 2. MOTIVATION
 3. BUG REPORTS and WISH LIST
 4. F2MATLAB CAPABILITIES
 5. F2MATLAB LIMITATIONS
 6. HOW TO USE F2MATLAB
 7. EXAMPLES
 8. REVISION HISTORY

-1.SUPPORT f2matlab. 
   I now also do conversion/translation/validation/optimization consulting.
   Please refer to my webpage:
   http://engineering.dartmouth.edu/~d30574x/consulting/consultingIndex.html

   Even though f2matlab is free (under GPL) for the using, I would like
   to ask that those who find it useful, wish to support the project,
   and are able to make a contribution to please do so commesurate with
   use (especially corporations). *** Important - Please donate using
   your PayPal account and not a credit card so as to avoid fees at
   PayPal. Thank you! PayPal email ID: barrowes@users.sourceforge.net


0. DISCLAIMER: Matlab is a trademark of the Mathworks company and is
   owned by them. The author makes no guarantee express or implied of
   any kind as to the applicability, usefulness, efficacy,
   bug-freeness, or the accuracy of the ensuing results from using
   f2matlab.

   The author bears no responsibility for any unwanted effect
   resulting from the use of this program. The author is not
   affiliated with the Mathworks.  The source code is given in full in
   the hopes that it will prove useful. Devlopment is done through
   sourceforge at f2matlab.sourceforge.net.

1. OBJECTIVE: f2matlab.m is a small translator which aims to
   convert Fortran90 files to Matlab m-files.

2. MOTIVATION: 

   1) Matlab is becoming ubiquitous in the engineering and scientific
   communities for its ease of use coupled with its powerful
   libraries. Yet the fact remains that a large number of stable and
   dependable programs exist in the fortran77/90 corpus. 

   2) Many times, often amidst the porting of fortran programs to
   Matlab, an automated converter of fortran90 code to Matlab code would
   be useful. 

   3) Having written matlab2fmex.m, a matlab to fortran90 mex file
   converter, the writing of f2matlab, which performs the reverse
   conversion, was substantially simplified.

3. BUG REPORTS and WISH LIST:
   For all bug reports, a wish list for f2matlab, and suggestions,
   see http://f2matlab.sourceforge.net/
   or email barrowes@users.sourceforge.net

4. F2MATLAB CAPABILITIES: f2matlab is aimed at converting 
   Fortran90 code to Matlab m-files. Accordingly, only basic data types
   and constructions are recommended. f2matlab can handle:

   all numeric types (handled by Matlab interpreter)
   most string functions
   comparisons, branches, loops, etc.
   basic read/write/print statements (if it's not too fancy...)
   modules

5. F2MATLAB LIMITATIONS: f2matlab can not handle some features of
   fortran90 yet. These include:

   can't handle complex read and write statements
   derived-typed variables
   equivalence
   ...

6. HOW TO USE F2MATLAB: f2matlab expects a single fortran90 fortran file to
   convert. If you have fortran77 code, use some free converter
   (e.g. to_f90 by Alan Miller) before running f2matlab. Then simply
   call f2matlab by using the full filename:
   f2matlab('filename.f90');
   The output will be filename.m in the same directory.
   
   A few flags are available that effect conversion:
   %  want_kb=0; 1 ==> if keyboard mode is desired after some conversion steps
   %  want_ze=0; 1 ==> direct f2matlab to zero all array variables.
   %  want_fi=0; 1 ==> direct f2matlab to try to put fix()'s around declared integers.

   Multiple subroutines and functions can and should be in the same fortran90 file.

7. EXAMPLES:
   Below are some specific examples. Other tests and examples include:
   
   Running the script TESTING_cfs in the comp_spec_func directory converts all of the f90
   code from Shanjie Zhang and Jianming Jin's book. These converted m-files can also be found
   at: http://ceta.mit.edu/comp_spec_func/

   For some larger examples (>=10000 lines), see the examples directory. The following
   fortran packages were downloaded (most from John Burkardt at FSU) and converted:

   - fempack routines for finite element analysis:
     f2matlab('femprb.f90');

   - quadrule
     f2matlab('quadrule_prb.f90');

   - quadrature integration routines by Arthur Stroud:
     f2matlab('stroud_prb.f90');


     Run the script TESTING_ex in the examples directory to convert and test all the examples
     there.



   The following f77 program, lpmns.f, calculates the associated legendre
   function for degree m, order n, and parameter x. The output is the function value pm
   and the value of the derivative pd. This function is from the book "Computation of
   Special Functions." by Shanjie Zhang, and Jianming Jin, New York : Wiley, c1996.

   These examples are in the examples directory.

       SUBROUTINE LPMNS(M,N,X,PM,PD)
c**************************************************************
C
C       ========================================================
C       Purpose: Compute associated Legendre functions Pmn(x)
C                and Pmn'(x) for a given order
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(n) --- Pmn(x)
C                PD(n) --- Pmn'(x)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PM(0:N),PD(0:N)
        DO 10 K=0,N
           PM(K)=0.0D0
10         PD(K)=0.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 K=0,N
              IF (M.EQ.0) THEN
                 PM(K)=1.0D0
                 PD(K)=0.5D0*K*(K+1.0)
                 IF (X.LT.0.0) THEN
                    PM(K)=(-1)**K*PM(K)
                    PD(K)=(-1)**(K+1)*PD(K)
                 ENDIF
              ELSE IF (M.EQ.1) THEN
                 PD(K)=1.0D+300
              ELSE IF (M.EQ.2) THEN
                 PD(K)=-0.25D0*(K+2.0)*(K+1.0)*K*(K-1.0)
                 IF (X.LT.0.0) PD(K)=(-1)**(K+1)*PD(K)
              ENDIF
15         CONTINUE
           RETURN
        ENDIF
        X0=DABS(1.0D0-X*X)
        PM0=1.0D0
        PMK=PM0
        DO 20 K=1,M
           PMK=(2.0D0*K-1.0D0)*DSQRT(X0)*PM0
20         PM0=PMK
        PM1=(2.0D0*M+1.0D0)*X*PM0
        PM(M)=PMK
        PM(M+1)=PM1
        DO 25 K=M+2,N
           PM2=((2.0D0*K-1.0D0)*X*PM1-(K+M-1.0D0)*PMK)/(K-M)
           PM(K)=PM2
           PMK=PM1
25         PM1=PM2
        PD(0)=((1.0D0-M)*PM(1)-X*PM(0))/(X*X-1.0)  
        DO 30 K=1,N
30          PD(K)=(K*X*PM(K)-(K+M)*PM(K-1))/(X*X-1.0D0)
        RETURN
        END

As this is fortran77, we use to_f90 to get lpmns.f90:

SUBROUTINE lpmns(m,n,x,pm,pd)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-14  Time: 13:32:33
 
!**************************************************************

!       ========================================================
!       Purpose: Compute associated Legendre functions Pmn(x)
!                and Pmn'(x) for a given order
!       Input :  x --- Argument of Pmn(x)
!                m --- Order of Pmn(x),  m = 0,1,2,...,n
!                n --- Degree of Pmn(x), n = 0,1,2,...,N
!       Output:  PM(n) --- Pmn(x)
!                PD(n) --- Pmn'(x)
!       ========================================================


INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(OUT)            :: pm(0:n)
DOUBLE PRECISION, INTENT(OUT)            :: pd(0:n)
IMPLICIT DOUBLE PRECISION (a-h,o-z)


DO  k=0,n
  pm(k)=0.0D0
  pd(k)=0.0D0
END DO
IF (DABS(x) == 1.0D0) THEN
  DO  k=0,n
    IF (m == 0) THEN
      pm(k)=1.0D0
      pd(k)=0.5D0*k*(k+1.0)
      IF (x < 0.0) THEN
        pm(k)=(-1)**k*pm(k)
        pd(k)=(-1)**(k+1)*pd(k)
      END IF
    ELSE IF (m == 1) THEN
      pd(k)=1.0D+300
    ELSE IF (m == 2) THEN
      pd(k)=-0.25D0*(k+2.0)*(k+1.0)*k*(k-1.0)
      IF (x < 0.0) pd(k)=(-1)**(k+1)*pd(k)
    END IF
  END DO
  RETURN
END IF
x0=DABS(1.0D0-x*x)
pm0=1.0D0
pmk=pm0
DO  k=1,m
  pmk=(2.0D0*k-1.0D0)*DSQRT(x0)*pm0
  pm0=pmk
END DO
pm1=(2.0D0*m+1.0D0)*x*pm0
pm(m)=pmk
pm(m+1)=pm1
DO  k=m+2,n
  pm2=((2.0D0*k-1.0D0)*x*pm1-(k+m-1.0D0)*pmk)/(k-m)
  pm(k)=pm2
  pmk=pm1
  pm1=pm2
END DO
pd(0)=((1.0D0-m)*pm(1)-x*pm(0))/(x*x-1.0)
DO  k=1,n
  pd(k)=(k*x*pm(k)-(k+m)*pm(k-1))/(x*x-1.0D0)
END DO
RETURN
END SUBROUTINE lpmns

Now we convert this to a Matlab function lpmns.m:

>> f2matlab('lpmns.f90');
  subroutine lpmns(m,n,x,pm,pd)
    Number of lines:   70
      Writing file:  lpmns.m ... completed 
  Total time: 0.58333
>> 

Resulting in lpmns.m:

function [m,n,x,pm,pd]=lpmns(m,n,x,pm,pd);

%**************************************************************

%       ========================================================
%       Purpose: Compute associated Legendre functions Pmn(x)
%                and Pmn'(x) for a given order
%       Input :  x --- Argument of Pmn(x)
%                m --- Order of Pmn(x),  m = 0,1,2,...,n
%                n --- Degree of Pmn(x), n = 0,1,2,...,N
%       Output:  PM(n) --- Pmn(x)
%                PD(n) --- Pmn'(x)
%       ========================================================


for  k=0:n;
pm(k+1)=0.0d0;
pd(k+1)=0.0d0;
end;
if (abs(x) == 1.0d0);
for  k=0:n;
if (m == 0);
pm(k+1)=1.0d0;
pd(k+1)=0.5d0.*k.*(k+1.0);
if (x < 0.0);
pm(k+1)=(-1).^k.*pm(k+1);
pd(k+1)=(-1).^(k+1).*pd(k+1);
end;
elseif (m == 1);
pd(k+1)=1.0d+300;
elseif (m == 2);
pd(k+1)=-0.25d0.*(k+2.0).*(k+1.0).*k.*(k-1.0);
if (x < 0.0) pd(k+1)=(-1).^(k+1).*pd(k+1); end;
end;
end;
return;
end;
x0=abs(1.0d0-x.*x);
pm0=1.0d0;
pmk=pm0;
for  k=1:m;
pmk=(2.0d0.*k-1.0d0).*sqrt(x0).*pm0;
pm0=pmk;
end;
pm1=(2.0d0.*m+1.0d0).*x.*pm0;
pm(m+1)=pmk;
pm(m+1+1)=pm1;
for  k=m+2:n;
pm2=((2.0d0.*k-1.0d0).*x.*pm1-(k+m-1.0d0).*pmk)./(k-m);
pm(k+1)=pm2;
pmk=pm1;
pm1=pm2;
end;
pd(0+1)=((1.0d0-m).*pm(1+1)-x.*pm(0+1))./(x.*x-1.0);
for  k=1:n;
pd(k+1)=(k.*x.*pm(k+1)-(k+m).*pm(k-1+1))./(x.*x-1.0d0);
end;




Now we can run lpmns.m:

>> [a1,a2,a3,pm,pd]=lpmns(0,10,.5)
a1 =
     0
a2 =
    10
a3 =
                       0.5
pm =
  Columns 1 through 4 
                         1                       0.5                    -0.125                   -0.4375
  Columns 5 through 8 
                -0.2890625                0.08984375              0.3232421875             0.22314453125
  Columns 9 through 11 
        -0.073638916015625        -0.267898559570312        -0.188228607177734
pd =
  Columns 1 through 4 
                         0                         1                       1.5                     0.375
  Columns 5 through 8 
                   -1.5625                -2.2265625               -0.57421875              1.9755859375
  Columns 9 through 11 
             2.77294921875         0.723724365234375         -2.31712341308594
>> 

Incidentally, it is now possible to convert lpmns.m into a fortran90
mex file callable from Matlab using matlab2fmex.m.
(See http://matlab2fmex.sourceforge.net/ )
First, construct another input which is the same size as the
output. Here is the modified lpmns2.m

function [pm,pd]=lpmns(m,n,x,outsize)

%**************************************************************

%       ========================================================
%       Purpose: Compute associated Legendre functions Pmn(x)
%                and Pmn'(x) for a given order
%       Input :  x --- Argument of Pmn(x)
%                m --- Order of Pmn(x),  m = 0,1,2,...,n
%                n --- Degree of Pmn(x), n = 0,1,2,...,N
%       Output:  PM(n) --- Pmn(x)
%                PD(n) --- Pmn'(x)
%       ========================================================

for  k=0:n;
 pm(k+1)=0.0d0;
 pd(k+1)=0.0d0;
end;
if (abs(x) == 1.0d0);
 for  k=0:n;
  if (m == 0);
   pm(k+1)=1.0d0;
   pd(k+1)=0.5d0.*k.*(k+1.0);
   if (x < 0.0);
    pm(k+1)=(-1).^k.*pm(k+1);
    pd(k+1)=(-1).^(k+1).*pd(k+1);
   end;
  elseif (m == 1);
   pd(k+1)=1.0d+300;
  elseif (m == 2);
   pd(k+1)=-0.25d0.*(k+2.0).*(k+1.0).*k.*(k-1.0);
   if (x < 0.0)
    pd(k+1)=(-1).^(k+1).*pd(k+1);
   end;
  end;
 end;
 return;
end;
x0=abs(1.0d0-x.*x);
pm0=1.0d0;
pmk=pm0;
for  k=1:m;
 pmk=(2.0d0.*k-1.0d0).*sqrt(x0).*pm0;
 pm0=pmk;
end;
pm1=(2.0d0.*m+1.0d0).*x.*pm0;
pm(m+1)=pmk;
pm(m+1+1)=pm1;
for  k=m+2:n;
 pm2=((2.0d0.*k-1.0d0).*x.*pm1-(k+m-1.0d0).*pmk)./(k-m);
 pm(k+1)=pm2;
 pmk=pm1;
 pm1=pm2;
end;
pd(0+1)=((1.0d0-m).*pm(1+1)-x.*pm(0+1))./(x.*x-1.0);
for  k=1:n;
 pd(k+1)=(k.*x.*pm(k+1)-(k+m).*pm(k-1+1))./(x.*x-1.0d0);
end;

Now create a lpmns.mat file using matlab2fmex_save:

>> n=10;matlab2fmex_save('[x,w]=lpmns2(0,n,.5,zeros(1,n+1));');

And now convert lpmns.m to a fortran90 mex file:

>> matlab2fmex('lpmns2');
Converting --- lpmns2.m  ==>  lpmns2.f  ==>  lpmns2.mex
matlab2fmex.
 ==> mex lpmns2.f mexfunctions.f mexoperators.f mexcallback.f 
>> 

Finally, we can call lpmns2 as a mex file from the Matlab command
prompt: 

>> [pm,pd]=lpmns2(0,n,.5,zeros(1,n+1))
pm =
  Columns 1 through 3 
                         1                       0.5                    -0.125
  Columns 4 through 6 
                   -0.4375                -0.2890625                0.08984375
  Columns 7 through 9 
              0.3232421875             0.22314453125        -0.073638916015625
  Columns 10 through 11 
        -0.267898559570312        -0.188228607177734
pd =
  Columns 1 through 3 
                         0                         1                       1.5
  Columns 4 through 6 
                     0.375                   -1.5625                -2.2265625
  Columns 7 through 9 
               -0.57421875              1.9755859375             2.77294921875
  Columns 10 through 11 
         0.723724365234375         -2.31712341308594
>> 
>> which lpmns2
/home/barrowes/compiler/f2matlab/examples/lpmns2.mexaxp
>> 

We have just converted a native fortran77 subroutine into a fortran90
mex file usable (with variable sized inputs) in Matlab.


lqmns.f in the examples directory can be handles similarly.


As another example, let's follow the same procedure with dgauleg.f, a
subroutine which calculates gaussian quadrature coefficients.

      SUBROUTINE dgauleg(x1,x2,x,w,n)
c----------------------------------------------------------------
      INTEGER n
      real*8 x1,x2,x(n),w(n)
      real*8 EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      real*8 p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      z1=huge(1.0)
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/
     & (n+.5d0))
        do while (abs(z-z1).gt.EPS)
         p1=1.d0
         p2=0.d0
         do 11 j=1,n
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 11      continue
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp
        end do
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end

Now convert to fortran90

SUBROUTINE dgauleg(x1,x2,x,w,n)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-14  Time: 13:40:52
 
!----------------------------------------------------------------

REAL*8, INTENT(IN OUT)                   :: x1
REAL*8, INTENT(IN OUT)                   :: x2
REAL*8, INTENT(OUT)                      :: x(n)
REAL*8, INTENT(OUT)                      :: w(n)
INTEGER, INTENT(IN)                      :: n


REAL*8
REAL*8, PARAMETER :: eps=3.d-14
INTEGER :: i,j,m
REAL*8 p1,p2,p3,pp,xl,xm,z,z1

m=(n+1)/2
z1=huge(1.0)
xm=0.5D0*(x2+x1)
xl=0.5D0*(x2-x1)
DO  i=1,m
  z=COS(3.141592654D0*(i-.25D0)/ (n+.5D0))
  DO WHILE (ABS(z-z1) > eps)
    p1=1.d0
    p2=0.d0
    DO  j=1,n
      p3=p2
      p2=p1
      p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
    END DO
    pp=n*(z*p1-p2)/(z*z-1.d0)
    z1=z
    z=z1-p1/pp
  END DO
  x(i)=xm-xl*z
  x(n+1-i)=xm+xl*z
  w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
  w(n+1-i)=w(i)
END DO
RETURN
END SUBROUTINE dgauleg

And now to dgauleg.m

>> f2matlab('dgauleg.f90');
  subroutine dgauleg(x1,x2,x,w,n)
    Number of lines:   44
  No Matlab function for 'huge' except replacing with 'realmax' in Matlab
    Writing file:  dgauleg.m ... completed 
  Total time: 0.35
>> 

Resulting dgauleg.m:

function [x1,x2,x,w,n]=dgauleg(x1,x2,x,w,n);

%----------------------------------------------------------------

eps=3.d-14;

m=(n+1)./2;
z1=realmax(1.0);
xm=0.5d0.*(x2+x1);
xl=0.5d0.*(x2-x1);
for  i=1:m;
z=cos(3.141592654d0.*(i-.25d0)./ (n+.5d0));
while (abs(z-z1) > eps);
p1=1.d0;
p2=0.d0;
for  j=1:n;
p3=p2;
p2=p1;
p1=((2.d0.*j-1.d0).*z.*p2-(j-1.d0).*p3)./j;
end;
pp=n.*(z.*p1-p2)./(z.*z-1.d0);
z1=z;
z=z1-p1./pp;
end;
x(i)=xm-xl.*z;
x(n+1-i)=xm+xl.*z;
w(i)=2.d0.*xl./((1.d0-z.*z).*pp.*pp);
w(n+1-i)=w(i);
end;


And now we can call it from Matlab:
>>  [x,w]=dgauleg(0,1,0,1,10)  
x =
     0
w =
     1
>>  [x1,x2,x,w,n]=dgauleg(0,1,0,1,10)
x1 =
     0
x2 =
     1
x =
  Columns 1 through 4 
        0.0130467357414142        0.0674683166555077         0.160295215850488         0.283302302935376
  Columns 5 through 8 
         0.425562830509184         0.574437169490816         0.716697697064624         0.839704784149512
  Columns 9 through 10 
         0.932531683344492         0.986953264258586
w =
  Columns 1 through 4 
        0.0333356721543414        0.0747256745752903         0.109543181257991         0.134633359654996
  Columns 5 through 8 
         0.147762112357376         0.147762112357376         0.134633359654996         0.109543181257991
  Columns 9 through 10 
        0.0747256745752903        0.0333356721543414
n =
    10
>> 


Finally, here is an example which converts a fortran90 program with subroutines into
Matlab. The program is MCPBDN.for from the book "Computation of Special Functions." by
Shanjie Zhang, and Jianming Jin, New York : Wiley, c1996. After using to_f90 on this file,
the only alteration was to specify the inputs x and y in the code instead of using a read
statement. See the /examples directory.

Converting:
>> f2matlab('mcpbdn.f90');
  program mcpbdn
    Number of lines:   61
  subroutine cpbdn(n,z,cpb,cpd)
    Number of lines:   103
  subroutine cpdsa(n,z,cdn)
    Number of lines:   57
  subroutine cpdla(n,z,cdn)
    Number of lines:   30
  subroutine gaih(x,ga)
    Number of lines:   29
    Writing file:  mcpbdn.m ... completed 
  Total time: 2.75
>> 

We can now try it out:

>> mcpbdn
ans =
please enter n, x and y 
ans =
  n     re[dn(z)]       im[dn(z)]       
ans =
re[dn'(z)]      im[dn'(z)]
ans =
-------------------------------------------
ans =
-------------------------
i =
     0
ans =
          0.997798279178581 +    0.0663218973512007i
ans =
          -2.32869095456845 -      2.66030044132445i
i =
     1
ans =
            4.6573819091369 +      5.32060088264891i
ans =
            2.6558457129586 -      24.8786350821133i
i =
     2
ans =
          -4.31389314673862 +      49.8235920615778i
ans =
           144.658476839065 -        103.1330455218i
i =
     3
ans =
          -280.002189859856 +      216.907292808898i
ans =
           1229.33202723167 +      307.208018812128i
i =
     4
ans =
          -2471.60573390356 -      464.945261439523i
ans =
           3896.64242172066 +      8209.00665959329i
i =
     5
ans =
          -8913.29360288074 -       15550.384147951i
ans =
          -28950.7550321934 +      58834.4680698817i
>> 

Which agrees with the values produced by the original code.



8. REVISION HISTORY:
   See changelog
