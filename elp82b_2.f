      subroutine ELP82B_2 (tjj,prec,nul,r,ierr)
*-----------------------------------------------------------------------
*
*     Reference : Bureau des Longitudes - MCJCGF.9601.
*
*     Object :
*     Computation of geocentric lunar coordinates from ELP 2000-82 and
*     ELP2000-85 theories (M. Chapront-Touze and J. Chapront).
*     Constants fitted to JPL's ephemerides DE200/LE200.
*
*     Remark :
*     This subroutine ELP82B_2 is used for computer memory larger than
*     1.3 Mo : Series are stocked in memory at the first call.
*     If the computer memory is too small, the subroutine ELP82B_1 has
*     to be used instead.
*
*     Input :
*     tjj    julian date TDB (real double precision).
*     prec   truncation level in radian (real double precision).
*     nul    number of logical unit for reading the files (integer).
*
*     Output :
*     r(3)   table of rectangular coordinates (real double precision).
*            reference frame : mean dynamical ecliptic and inertial
*            equinox of J2000 (JD 2451545.0).
*            r(1) : X (kilometer).
*            r(2) : Y (kilometer).
*            r(3) : Z (kilometer).
*     ierr   error index (integer).
*            ierr=0 : no error.
*            ierr=1 : error in elp 2000-82 files (end of file).
*            ierr=2 : error in elp 2000-82 files (reading error).
*
*     Files :
*     36 data files include the series related to various components of
*     the theory for the 3 spherical coordinates : longitude, latitude
*     and distance.
*
*     Documentation :
*     Files, series, constants and coordinate systems are described in
*     the notice LUNAR SOLUTION ELP 2000-82B.
*
*-----------------------------------------------------------------------
*
*     Declarations and initializations
*     --------------------------------
*
      implicit double precision (a-h,o-z)
      character fich*60
*
      dimension w(3,0:4),eart(0:4),peri(0:4),p(8,0:1)
      dimension del(4,0:4),zeta(0:1),t(0:4)
      dimension r(3),pre(3),coef(7),ilu(4),ipla(11)
      dimension nterm(3,12),nrang(3,12),zone(6)
      dimension pc1(6,1023),pc2(6,918),pc3(6,704)
      dimension per1(3,19537),per2(3,6766),per3(3,8924)
*
      parameter (cpi=3.141592653589793d0,cpi2=2.d0*cpi,pis2=cpi/2.d0)
      parameter (rad=648000.d0/cpi,deg=cpi/180.d0,c1=60.d0,c2=3600.d0)
      parameter (ath=384747.9806743165d0,a0=384747.9806448954d0)
      parameter (dj2000=2451545.0d0,sc=36525.d0)
*
      data ideb/0/,prec0/-1.d0/,t/1.d0,4*0.d0/
*
      r(1)=0.d0
      r(2)=0.d0
      r(3)=0.d0
      t(1)=(tjj-dj2000)/sc
      t(2)=t(1)*t(1)
      t(3)=t(2)*t(1)
      t(4)=t(3)*t(1)
*
*     Parameters
*     ----------
*
      if (ideb.eq.0) then
*
         ideb=1
*
         am=0.074801329518d0
         alpha=0.002571881335d0
         dtasm=2.d0*alpha/(3.d0*am)
*
*        Lunar arguments.
*
         w(1,0)=(218+18/c1+59.95571d0/c2)*deg
         w(2,0)=(83+21/c1+11.67475d0/c2)*deg
         w(3,0)=(125+2/c1+40.39816d0/c2)*deg
         eart(0)=(100+27/c1+59.22059d0/c2)*deg
         peri(0)=(102+56/c1+14.42753d0/c2)*deg
         w(1,1)=1732559343.73604d0/rad
         w(2,1)=14643420.2632d0/rad
         w(3,1)=-6967919.3622d0/rad
         eart(1)=129597742.2758d0/rad
         peri(1)=1161.2283d0/rad
         w(1,2)=-5.8883d0/rad
         w(2,2)=-38.2776d0/rad
         w(3,2)=6.3622d0/rad
         eart(2)=-0.0202d0/rad
         peri(2)=0.5327d0/rad
         w(1,3)=0.6604d-2/rad
         w(2,3)=-0.45047d-1/rad
         w(3,3)=0.7625d-2/rad
         eart(3)=0.9d-5/rad
         peri(3)=-0.138d-3/rad
         w(1,4)=-0.3169d-4/rad
         w(2,4)=0.21301d-3/rad
         w(3,4)=-0.3586d-4/rad
         eart(4)=0.15d-6/rad
         peri(4)=0.d0
*
*        Precession constant.
*
         precess=5029.0966d0/rad
*
*        Planetary arguments.
*
         p(1,0)=(252+15/c1+3.25986d0/c2)*deg
         p(2,0)=(181+58/c1+47.28305d0/c2)*deg
         p(3,0)=eart(0)
         p(4,0)=(355+25/c1+59.78866d0/c2)*deg
         p(5,0)=(34+21/c1+5.34212d0/c2)*deg
         p(6,0)=(50+4/c1+38.89694d0/c2)*deg
         p(7,0)=(314+3/c1+18.01841d0/c2)*deg
         p(8,0)=(304+20/c1+55.19575d0/c2)*deg
         p(1,1)=538101628.68898d0/rad
         p(2,1)=210664136.43355d0/rad
         p(3,1)=eart(1)
         p(4,1)=68905077.59284d0/rad
         p(5,1)=10925660.42861d0/rad
         p(6,1)=4399609.65932d0/rad
         p(7,1)=1542481.19393d0/rad
         p(8,1)=786550.32074d0/rad
*
*        Corrections of the constants (fit to DE200/LE200).
*
         delnu=+0.55604d0/rad/w(1,1)
         dele=+0.01789d0/rad
         delg=-0.08066d0/rad
         delnp=-0.06424d0/rad/w(1,1)
         delep=-0.12879d0/rad
*
*        Delaunay's arguments.
*
         do i=0,4
            del(1,i)=w(1,i)-eart(i)
            del(4,i)=w(1,i)-w(3,i)
            del(3,i)=w(1,i)-w(2,i)
            del(2,i)=eart(i)-peri(i)
         enddo
         del(1,0)=del(1,0)+cpi
         zeta(0)=w(1,0)
         zeta(1)=w(1,1)+precess
*
*        Precession matrix.
*
         p1=0.10180391d-4
         p2=0.47020439d-6
         p3=-0.5417367d-9
         p4=-0.2507948d-11
         p5=0.463486d-14
         q1=-0.113469002d-3
         q2=0.12372674d-6
         q3=0.1265417d-8
         q4=-0.1371808d-11
         q5=-0.320334d-14
*
      endif
*
*     Reading files
*     -------------
*
      if (prec.ne.prec0) then
*
         prec0=prec
         pre(1)=prec*rad-1.d-12
         pre(2)=prec*rad-1.d-12
         pre(3)=prec*ath
*
         do ific=1,36
*
            ir=0
            itab=(ific+2)/3
            iv=mod(ific-1,3)+1
*
            select case (ific)
*
*           Files : Main problem.
*
            case (1:3)
*
            write (fich,2000) ific
            open (nul,file=fich,err=600)
            read (nul,1000,end=500,err=600)
*
100         continue
            read (nul,1001,end=400,err=600) ilu,coef
            xx=coef(1)
            if (abs(xx).lt.pre(iv)) goto 100
*
            ir=ir+1
            tgv=coef(2)+dtasm*coef(6)
            if (ific.eq.3) coef(1)=coef(1)-2.d0*coef(1)*delnu/3.d0
            xx=coef(1)+tgv*(delnp-am*delnu)+coef(3)*delg+coef(4)*dele
     .        +coef(5)*delep
            zone(1)=xx
            do k=0,4
               y=0.d0
               do i=1,4
                  y=y+ilu(i)*del(i,k)
              enddo
               zone(k+2)=y
            enddo
            if (iv.eq.3) zone(2)=zone(2)+pis2
            do i=1,6
               if (iv.eq.1) pc1(i,ir)=zone(i)
               if (iv.eq.2) pc2(i,ir)=zone(i)
               if (iv.eq.3) pc3(i,ir)=zone(i)
            enddo
*
            goto 100
*
*           Files : Tides - Relativity - Solar eccentricity.
*
            case (4:9,22:36)
*
            if (ific.le.9) write (fich,2000) ific
            if (ific.gt.9) write (fich,2001) ific
            open (nul,file=fich,err=600)
            read (nul,1000,end=500,err=600)
*
200         continue
            read (nul,1002,end=400,err=600) iz,ilu,pha,xx
            if (xx.lt.pre(iv)) goto 200
*
            ir=ir+1
            zone(1)=xx
            do k=0,1
               if (k.eq.0) y=pha*deg
               if (k.ne.0) y=0.d0
               y=y+iz*zeta(k)
               do i=1,4
                  y=y+ilu(i)*del(i,k)
              enddo
              zone(k+2)=y
            enddo
            j=nrang(iv,itab-1)+ir
            do i=1,3
               if (iv.eq.1) per1(i,j)=zone(i)
               if (iv.eq.2) per2(i,j)=zone(i)
               if (iv.eq.3) per3(i,j)=zone(i)
            enddo
*
            goto 200
*
*           Files : Planetary perturbations.
*
            case (10:21)
*
            write (fich,2001) ific
            open (nul,file=fich,err=600)
            read (nul,1000,end=500,err=600)
*
300         continue
            read (nul,1003,end=400,err=600) ipla,pha,xx
            if (xx.lt.pre(iv)) goto 300
*
            ir=ir+1
            zone(1)=xx
            if (ific.lt.16) then
               do k=0,1
                  if (k.eq.0) y=pha*deg
                  if (k.ne.0) y=0.d0
                  y=y+ipla(9)*del(1,k)+ipla(10)*del(3,k)
     .               +ipla(11)*del(4,k)
                  do i=1,8
                     y=y+ipla(i)*p(i,k)
                  enddo
                  zone(k+2)=y
               enddo
            else
               do k=0,1
                  if (k.eq.0) y=pha*deg
                  if (k.ne.0) y=0.d0
                  do i=1,4
                     y=y+ipla(i+7)*del(i,k)
                  enddo
                  do i=1,7
                     y=y+ipla(i)*p(i,k)
                  enddo
                  zone(k+2)=y
               enddo
            endif
            j=nrang(iv,itab-1)+ir
            do i=1,3
               if (iv.eq.1) per1(i,j)=zone(i)
               if (iv.eq.2) per2(i,j)=zone(i)
               if (iv.eq.3) per3(i,j)=zone(i)
            enddo
*
            goto 300
*
            endselect
*
400         continue
*
            nterm(iv,itab)=ir
            if (itab.eq.1) then
               nrang(iv,itab)=0
            else
               nrang(iv,itab)=nrang(iv,itab-1)+nterm(iv,itab)
            endif
            close (nul)
*
         enddo
*
      endif
*
*     Substitution of time
*     --------------------
*
      do iv=1,3
         r(iv)=0.d0
*
         do itab=1,12
*
            do nt=1,nterm(iv,itab)
*
               select case (iv)
*
                  case (1)
                  if (itab.eq.1) then
                     x=pc1(1,nt)
                     y=pc1(2,nt)
                     do k=1,4
                        y=y+pc1(k+2,nt)*t(k)
                     enddo
                  else
                     j=nrang(1,itab-1)+nt
                     x=per1(1,j)
                     y=per1(2,j)+per1(3,j)*t(1)
                  endif
*
                  case (2)
                  if (itab.eq.1) then
                     x=pc2(1,nt)
                     y=pc2(2,nt)
                     do k=1,4
                        y=y+pc2(k+2,nt)*t(k)
                     enddo
                  else
                     j=nrang(2,itab-1)+nt
                     x=per2(1,j)
                     y=per2(2,j)+per2(3,j)*t(1)
                  endif
*
                  case (3)
                  if (itab.eq.1) then
                     x=pc3(1,nt)
                     y=pc3(2,nt)
                     do k=1,4
                        y=y+pc3(k+2,nt)*t(k)
                     enddo
                  else
                     j=nrang(3,itab-1)+nt
                     x=per3(1,j)
                     y=per3(2,j)+per3(3,j)*t(1)
                  endif
*
               endselect
*
               if (itab.eq.3.or.itab.eq.5.or.
     .             itab.eq.7.or.itab.eq.9) x=x*t(1)
               if (itab.eq.12) x=x*t(2)
*
               r(iv)=r(iv)+x*sin(y)
*
            enddo
*
         enddo
*
      enddo
*
*     Change of coordinates
*     ---------------------
*
      r(1)=r(1)/rad+w(1,0)+w(1,1)*t(1)+w(1,2)*t(2)+w(1,3)*t(3)
     .    +w(1,4)*t(4)
      r(2)=r(2)/rad
      r(3)=r(3)*a0/ath
*
      x1=r(3)*cos(r(2))
      x2=x1*sin(r(1))
      x1=x1*cos(r(1))
      x3=r(3)*sin(r(2))
*
      pw=(p1+p2*t(1)+p3*t(2)+p4*t(3)+p5*t(4))*t(1)
      qw=(q1+q2*t(1)+q3*t(2)+q4*t(3)+q5*t(4))*t(1)
      ra=2.d0*sqrt(1.d0-pw*pw-qw*qw)
      pwqw=2.d0*pw*qw
      pw2=1.d0-2*pw*pw
      qw2=1.d0-2*qw*qw
      pw=pw*ra
      qw=qw*ra
*
      r(1)=pw2*x1+pwqw*x2+pw*x3
      r(2)=pwqw*x1+qw2*x2-qw*x3
      r(3)=-pw*x1+qw*x2+(pw2+qw2-1.d0)*x3
*
      ierr=0
      return
*
*
*     Errors
*     ------
*
500   continue
*
      ierr=1
      return
*
600   continue
*
      ierr=2
      return
*
*     Formats
*     -------
*
1000  format (1x)
1001  format (4i3,2x,f13.5,6(2x,f10.2))
1002  format (5i3,1x,f9.5,1x,f9.5)
1003  format (11i3,1x,f9.5,1x,f9.5)
2000  format ('ELP',i1,1x)
2001  format ('ELP',i2)
*
        end
