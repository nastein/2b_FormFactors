program test_iso

   use dirac_matrices
   implicit none

   integer*4 :: i,j,k,l
   real*8 :: test_p4(4),w 
   complex*16 :: t1(2), t1p(2), t2(2), t2p(2), p(2), n(2), result, result2
   complex*16 :: czero = (0.0d0,0.0d0)
   complex*16 :: cone  = (1.0d0,0.0d0)
   complex*16 :: ci    = (0.0d0,1.0d0)

   w = 0.0d0
   test_p4 = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)


   call dirac_matrices_in(938.0d0,938.0d0,100.0d0,0.0d0,0.0d0)
   call current_init(w,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,2)
   call define_spinors()

   p(1) = cone
   p(2) = czero

   n(1) = czero
   n(2) = cone

   t1 = n
   t2 = p
   t1p = p 
   t2p = p

   !Ivminus is called Ivminus(it1,it2,itp1,itp2)

   result = -Ivminus(p,p,n,p)
   !write(6,*) '<n1p2 | Iv^dag | p2p1> = ', result

   result = -Ivminus(n,p,n,n)
   !write(6,*) '<n1n2 | Iv^dag | n2p1> = ', result

   Write(6,*) '12 -> 21 exchange:'
      result2 = IDeltaBdag(p,p,n,p)
   write(6,*) '<n1p2 | DeltaB^dag | p2p1> = ', result2

   result2 = IDeltaBdag(n,p,n,n)
   write(6,*) '<n1n2 | DeltaB^dag | n2p1> = ', result2

   Write(6,*) '21 -> 12 exchange:'
   result2 = IDeltaDdag(p,p,p,n)
   write(6,*) '<p2n1 | DeltaD^dag | p1p2> = ', result2

   result2 = IDeltaDdag(p,n,n,n)
   write(6,*) '<n2n1 | DeltaD^dag | p1n2> = ', result2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Write(6,*) '12 -> 21 exchange:'
   result2 = IDeltaDdag(p,p,n,p)
   write(6,*) '<n1p2 | DeltaD^dag | p2p1> = ', result2

   result2 = IDeltaDdag(n,p,n,n)
   write(6,*) '<n1n2 | DeltaD^dag | n2p1> = ', result2

   Write(6,*) '21 -> 12 exchange:'
   result2 = IDeltaBdag(p,p,p,n)
   write(6,*) '<p2n1 | DeltaB^dag | p1p2> = ', result2

   result2 = IDeltaBdag(p,n,n,n)
   write(6,*) '<n2n1 | DeltaB^dag | p1n2> = ', result2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Write(6,*) '12 -> 21 exchange:'
   result2 = IDeltaCdag(p,p,n,p)
   write(6,*) '<n1p2 | DeltaC^dag | p2p1> = ', result2

   result2 = IDeltaCdag(n,p,n,n)
   write(6,*) '<n1n2 | DeltaC^dag | n2p1> = ', result2

   Write(6,*) '21 -> 12 exchange:'
   result2 = IDeltaAdag(p,p,p,n)
   write(6,*) '<p2n1 | DeltaA^dag | p1p2> = ', result2

   result2 = IDeltaAdag(p,n,n,n)
   write(6,*) '<n2n1 | DeltaA^dag | p1n2> = ', result2


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Write(6,*) '12 -> 21 exchange:'
   result2 = IDeltaAdag(p,p,n,p)
   write(6,*) '<n1p2 | DeltaA^dag | p2p1> = ', result2

   result2 = IDeltaAdag(n,p,n,n)
   write(6,*) '<n1n2 | DeltaA^dag | n2p1> = ', result2

   Write(6,*) '21 -> 12 exchange:'
   result2 = IDeltaCdag(p,p,p,n)
   write(6,*) '<p2n1 | DeltaC^dag | p1p2> = ', result2

   result2 = IDeltaCdag(p,n,n,n)
   write(6,*) '<n2n1 | DeltaC^dag | p1n2> = ', result2



   

   


   

   
end program test_iso


