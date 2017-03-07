program md  !simple MD program
  implicit none
  double precision:: deltat, t, tmax, vx(500), vy(500), xcord(500), ycord(500),deltax(500), deltay(500),f(500)
  double precision::  deltar(500),kenergy(500), totkenergy
  double precision,parameter:: mass=39.948, sigma=3.405, epsillon=0.01035, npart=500, lx=60 ,ly=60   
  real,external :: ran3

  write(*,*) '################################ welcome to MD simulation ###############################'
  write(*,*) '############################### 500 Argon particle in in cell #########################'

  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) '###############################    simulation start    #########################################'
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''

  
  call init(xcord, ycord, vx, vy,deltax, deltay, deltar,totkenergy)  !initialization
  
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''

  
   !t=0
   !deltat=0
   !tmax=0.00000000001
   !write(*,*) tmax
   !deltat=tmax/2000
   !write(*,*) deltat

  t=0
  deltat=0
  tmax=20
  deltat=2

  
  do while (t<tmax) !MD loop
     call force(xcord, ycord, vx, vy,deltax, deltay, deltar,totkenergy,f)  !determine the forces
     !   call integrate(f,en) !integrate equations of motion
     t=t+deltat
     !write(*,*) t
     !  call sample  !sample averages
  enddo
  !  stop
  
  
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) '#####################################  simulation finish  ######################################'
end program md


!#############################################################################################################
!#############################################################################################################


subroutine init(xcord, ycord, vx, vy,deltax, deltay, deltar,totkenergy)
  implicit none
  double precision:: deltat, t, tmax, vx(500), vy(500), xcord(500), ycord(500),deltax(500), deltay(500), f(500)
  double precision::  deltar(500),kenergy(500), totkenergy
  double precision,parameter:: mass=39.948, sigma=3.405, epsillon=0.01035, npart=500,lx=60,ly=60   
  real,external :: ran3
  INTEGER:: idum,i
  real:: r1,u

  write(*,*) '#################################### The Initialize Start ######################################'

  xcord=0
  ycord=0
  deltax=0
  deltay=0
  deltar=0
  vx=0
  vy=0
  kenergy=0
  totkenergy=0
  
  OPEN(1,file='xcord.txt')
  OPEN(2,file='ycord.txt')
  OPEN(3,file='vx.txt')
  OPEN(4,file='vy.txt')
  OPEN(5,file='total.txt')
!  OPEN(6,file='kenergy.txt')
  
  
  do  i=1,npart
     
     idum= i*100+2000
     xcord(i)=ran3(idum)*60 !this scale to 0 to 60
     write(1,*) i, xcord(i)
     write(5,*) i, xcord(i)
     
     idum= i*100+1000
     ycord(i)=ran3(idum)*60  !this scale to 0 to 60
     write(2,*) i, ycord(i)
     write(5,*) ycord(i)
     
     idum=i*100+500
     vx(i)=ran3(idum)
     write(3,*) i,vx(i)
     write(5,*)  vx(i)
     
     idum=i*150+600
     vy(i)=ran3(idum)
     write(4,*) i,vy(i)
     write(5,*) vy(i)
     
     kenergy(i)=0.5*mass*(vx(i)**2+vy(i)**2)
     ! write(6,*) i, kenergy(i)
     
     totkenergy=totkenergy+kenergy(i)
  enddo
  
  !  totkenergy=totkenergy/npart
  
  write(*,*) ''
  write(*,*) '' 
  write(*,*) 'total kinetic energy:'
  write(*,*) totkenergy

  write(*,*) ''
  write(*,*) ''

  write(*,*) '############################################# The Initialize Finish #############################'
  
end subroutine init

!#############################################################################################################
!#############################################################################################################



subroutine force(xcord, ycord, vx, vy,deltax, deltay, deltar,totkenergy,f)
  implicit none
  double precision:: deltat, t, tmax, vx(500), vy(500), xcord(500), ycord(500),deltax(500), deltay(500), f(500)
  double precision::  deltar(500),kenergy(500), totkenergy
  double precision,parameter:: mass=39.948, sigma=3.405, epsillon=0.01035, npart=500, lx=60 ,ly=60 
  INTEGER:: i, j , k
  double precision:: rm2,rm6,rm12,phi,dphi
  
  do k=1,npart
     f(k)=0            !set force to zero
  enddo

  OPEN(7,file='deltaghablperiod.txt')
  OPEN(8,file='deltabadeperiod.txt')
  !OPEN(5,file='total.txt')
  
  i=0
  do i=1,npart-1
     write(*,*)'########################################### the i chanded ###########################################'
     write(*,*)'########################################### the i chanded ###########################################'
     write(*,*)'########################################### the i chanded ###########################################'
     write(*,*)'########################################### the i chanded ###########################################'
     write(*,*)'########################################### the i chanded ###########################################'
     do j=i+1, npart
        write(*,*)'###############   the j equal   ##################'
        write (*,*) 'ghabl'
        deltax(j)=xcord(i)-xcord(j)
        deltay(j)=ycord(i)-ycord(j)        
        ! write (7,*) i, j, deltax(j)
        write (*,*) deltax(j)
        
        write (*,*)'deltay'
        write (*,*) deltay(j)
        !periodic condition in x 
        
        if (deltax(j) < -30) then
           deltax(j)= deltax(j) +60
           write(*,*) ''
           write(*,*) ''
           write(*,*) 'manfi'
           write(*,*) ''
           write(*,*) ''
           write (*,*) deltax(j)
        elseif(deltax(j) > 30) then
           deltax(j)= deltax(j)-60
           write(*,*) ''
           write(*,*) ''
           write(*,*) 'mosbat'
           write(*,*) ''
           write(*,*) ''
           write (*,*) deltax(j)
        else
           write(*,*) ''
           write(*,*) ''
           write(*,*) 'bedone taghir'
           write(*,*) ''
           write(*,*) ''
           write (*,*) deltax(j)   
        endif
        !periodic condition in y

        write (*,*)'deltay'
        
        if (deltay(j) < -30) then
           deltay(j)= deltay(j) +60
           write(*,*) ''
           write(*,*) ''
           write(*,*) 'manfi'
           write(*,*) ''
           write(*,*) ''
           write (*,*) deltay(j)
        elseif(deltay(j) > 30) then
           deltay(j)= deltay(j)-60
           write(*,*) ''
           write(*,*) ''
           write(*,*) 'mosbat'
           write(*,*) ''
           write(*,*) ''
           write(*,*) deltay(j)
        else
           write(*,*) ''
           write(*,*) ''
           write(*,*) 'bedone taghir'
           write(*,*) ''
           write(*,*) ''
           write(*,*) deltay(j)   
        endif
      !  write (8,*) j,deltax(j)
        !  deltay=ycord(i)-ycord(j)
        deltar(j)= sqrt(deltax(j)**2+deltay(j)**2)
        
        write(*,*) ''
        write(*,*) ''
        write(*,*) '2.5 sigma'
        write(*,*) ''
        write(*,*) ''
        write(*,*) 2.5*sigma 
        


        write(*,*) ''
        write(*,*) ''
        write(*,*) 'deltar'
        write(*,*) ''
        write(*,*) ''
        write(*,*) deltar(j) 
        
        if(deltar(j) <2.5*sigma) then   ! compute Lennard-Jones potntl
           write(*,*) "mohasebe shavad"
           rm2=1./deltar(j)  !1/r^2
           rm6=rm2**3        !1/r^6
           rm12=rm2**6       !1/r^12
           phi=4*epsillon*(rm12-rm6)             !4[1/r^12 - 1/r^6] - phi(Rc)
                                                 !  The following is dphi = -(1/r)(dV/dr)
           dphi=24.0*rm2*( 2.d0*rm12 - rm6 )     !  24[2/r^14 - 1/r^8]

           !print("hi")phi,
           write(*,*) dphi 
           
!C
!           rm2 = 1.d0/Rsqij                         !  1/r^2
!           rm6 = rm2**3                             !  1/r^6
!           rm12 = rm6**2                            !  1/r^12
!           phi  = 4.d0 * ( rm12 - rm6 ) - phicutoff !  4[1/r^12 - 1/r^6] - phi(Rc)
!                                      !  The following is dphi = -(1/r)(dV/dr)
!         dphi = 24.d0*rm2*( 2.d0*rm12 - rm6 )     !  24[2/r^14 - 1/r^8]
!         ene_pot(i) = ene_pot(i) + 0.5d0*phi      ! accumulate energy
!         ene_pot(j) = ene_pot(j) + 0.5d0*phi      ! (i and j share it)
!         virial = virial - dphi*Rsqij             ! accumul. virial=sum r(dV/dr)
!         acc(:,i) = acc(:,i) + dphi*Sij           ! accumulate forces
!         acc(:,j) = acc(:,j) - dphi*Sij           !    (Fji = -Fij)
!      endif
        else
            write(*,*) "mohasebe nashavad"
        end if
     enddo
  enddo
  ! if
  
end subroutine force

!#############################################################################################################
!#############################################################################################################     
!subroutine integrate
!end subroutine integrate
!#############################################################################################################
!#############################################################################################################
FUNCTION ran3(idum)
  INTEGER idum
  INTEGER MBIG,MSEED,MZ
  ! C     REAL MBIG,MSEED,MZ
  REAL ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
  ! C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
  INTEGER i,iff,ii,inext,inextp,k
  INTEGER mj,mk,ma(55)
  ! C     REAL mj,mk,ma(55)
  SAVE iff,inext,inextp,ma
  DATA iff /0/
  if(idum.lt.0.or.iff.eq.0)then
     iff=1
     mj=abs(MSEED-abs(idum))
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do 11 i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(ii)
11      continue
        do 13 k=1,4
           do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12            continue
13            continue
              inext=0
              inextp=31
              idum=1
           endif
           inext=inext+1
           if(inext.eq.56)inext=1
           inextp=inextp+1
           if(inextp.eq.56)inextp=1
           mj=ma(inext)-ma(inextp)
           if(mj.lt.MZ)mj=mj+MBIG
           ma(inext)=mj
           ran3=mj*FAC
           return
        
end FUNCTION ran3       
!#############################################################################################################
!#############################################################################################################
