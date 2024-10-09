program test_iso

   use mag_mom_current
   implicit none
   integer*4 :: i,j
   real*8 :: k(4),kp(4),q(4),kp_3(3),kp_mag,mN,E
   complex*16 :: czero = (0.0d0,0.0d0)
   complex*16 :: cone  = (1.0d0,0.0d0)
   complex*16 :: ci    = (0.0d0,1.0d0)
   complex*16 :: tens(4,4)

   k = (/2.0d0, 0.0d0, 0.0d0, 2.0d0/)
   kp_3 = (/0.0d0, 0.5d0, 0.3d0/)
   kp_mag = sqrt(kp_3(1)**2 + kp_3(2)**2 + kp_3(3)**2)
   mN = 0.5d0
   E = sqrt(kp_mag**2 + mN**2)
   kp(1) = E 
   kp(2:4) = kp_3(:)

   q = k - kp
   write(6,*)'k = ', k 
   write(6,*)'kp =', kp
   write(6,*)'q = ', q

   call dirac_matrices_in(0.0d0, mN)
   call lepton_current_init(k, kp, q)
   call define_lept_spinors()
   call curr(tens)

   do i=1,4
      do j=1,4
         write(6,*)tens(i,j)
      enddo
   enddo

   
   
end program test_iso


