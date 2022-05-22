      subroutine coul_gen_cc_rel(Z,ee,r,lop,doprint,f,g,fd,gd,sig)
      implicit none
      integer lop,i,j,lrange,ifail
      real*8 Z,r,fsc_inv
      complex*16 eta,xlmin,k,rk,sig(1)
      complex*16 f(1),g(1),fd(1),gd(1)
      complex*16 ee
      logical doprint
	  !Variable fsc_inv added to calculate relativistic k
      ifail=0
      lrange=1
      xlmin=dcmplx(dble(lop),0.d0)

	  fsc_inv = 137.035999074
      k=sqrt(2.d0*ee)*sqrt(1.d0+ee/(2.d0*fsc_inv**2.d0))    ! relativistic k
      eta=-Z*(1.d0+ee/fsc_inv**2.d0)/k                      ! relativistic eta
* with real eta Re(f)= regular sol Re(g)= irregular sol
* the Im(f)=Re(g)
      call COULcc(r*k, ETA,xlmin,LRANGE, F,G,FD,GD,sig,1,0,IFAIL)
      if(ifail.ne.0) write(6,*)' Warning in Coulcc IFAIL=',ifail
      if(doprint.or.ifail.ne.0) then
        write
     :  (6,'(a,f5.1,a,1p2d12.4,a,1pd12.4,a,1p2d16.8,3(a,1p2d16.8))')
     :   ' Regular Coulomb function for Z=',Z,' E=',ee,
     :  ' at r=',r,' =',f(1),' d/dr',fd(1),' irregular',
     :   g(1),' d/de',gd(1)
        write(6,'(a,1p2d12.4)') ' sig=',sig(1)
      end if
      if(ifail.ne.0) then
        write(6,'(a)')
     : ' The program stops since Coulcc failed in coul_gen_cc'
        stop
      end if


      return
      end
