      subroutine b2xvsg (n, fun, inc, vrname, chsg)
      use b2mod_types
      implicit none
      integer n, inc
      real (kind=R8) ::
     *   fun(0:*)
      character vrname*(*), chsg*(*)
*     -----------------------------------------------------------------
*     B2XVSG verifies the sign of some vector quantity.
*     The vector fun is represented by (n,fun,inc), using BLAS storage
*     conventions. The sign of each element of fun must be as specified
*     by chsg, which may be '.lt.', '.le.', '.gt.' or '.ge.'. If this
*     is violated the code will halt via a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with fun.
*     -----------------------------------------------------------------
!      integer, save :: ncall = 0, no_test_sign = 1                       !srv 26.10.08  11.01.13
      character vrnam0*8
!!!      real (kind=R8) :: sminmax
      real (kind=R8) ::
     *   b2smin, b2smax
      external xertst, xerrab, b2smin, b2smax
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvsg')
*   ..set internal parameters on first call
!      if (ncall.eq.0) then
!       call ipgeti ('b2xvsg_no_test_sign', no_test_sign)                 !srv 26.10.08
!	 write(*,*) ' ncall =', ncall, ' no_test_sign =',no_test_sign
!      endif
      call xertst (0.lt.n.and.inc.ne.0,
     &  'b2xvsg--faulty argument n, inc')
      vrnam0 = vrname
!	if(no_test_sign.eq.0) then                                         !srv 26.10.08
!!!      if (vrnam0 .eq. 'zn') then
!!!         write(*,*) 'DPC: debug -- b2xvsg ns = ', n
!!!         write(*,*) 'DPC: debug -- b2xvsg zn = ', fun(0:n-1:inc)
!!!         write(*,*) 'DPC: debug -- b2xvsg chsg = ', chsg
!!!         sminmax = b2smin(n,fun,inc)
!!!         write(*,*) 'DPC: debug -- b2xvsg b2smin(n,fun,inc) = ',
!!!     &        sminmax
!!!      endif
      if (chsg.eq.'.lt.') then
        call xertst (b2smax(n,fun,inc).lt.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.le.') then
        call xertst (b2smax(n,fun,inc).le.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.gt.') then
!!!        write(*,*) 'DPC: debug -- b2xvsg vrnam0 = ', vrnam0
!!!        write(*,*) 'DPC: debug -- b2xvsg chsg = ', chsg
!!!        sminmax = b2smin(n,fun,inc)
!!!        write(*,*) 'DPC: debug -- b2xvsg (b2smin(n,fun,inc).gt.0) = ',
!!!     &        (sminmax.gt.0)
        call xertst (b2smin(n,fun,inc).gt.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.ge.') then
        call xertst (b2smin(n,fun,inc).ge.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else
        call xerrab ('b2xvsg--faulty argument chsg')
      endif
!	endif
!	ncall = ncall + 1
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvsg
