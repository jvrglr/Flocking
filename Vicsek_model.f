      program Vicksek_model
c        atan2 documentation--> https://gcc.gnu.org/onlinedocs/gcc-4.2.3/gfortran/ATAN2.html

        implicit none
        integer*4 N,L,t
        double precision PI,mean,s_m,c_m,r,v
        parameter (N=3,L=25,PI=4.0d0*datan(1.0d0),r=1.0d0,v=3.0d-1) !Number of birds
        double precision angle(N) !orientation of each bird, angle in [-PI,PI]
        double precision x(N),y(N) !x,y position of each bird
        double precision a_noise,amp !angular noise and amplitude of noise
        double precision dran_u
        integer*4 i,aux
        call dran_ini(1994)
        open(unit=1, file="data.dat", status="unknown")
        t=0
        amp=1.0
        call initialize(N,L,PI,x,y,angle)
        do i=1,N
          write(1,*) x(i),y(i),cos(angle(i)),sin(angle(i))
        enddo
        aux=2
        s_m=0.0d0
        c_m=0.0d0
        do i=1,N
          s_m=s_m+sin(angle(i))
          c_m=c_m+cos(angle(i))
        enddo
        s_m=s_m/dble(N)
        c_m=c_m/dble(N)
        angle=atan2(s_m,c_m)
        do i=1,N
          angle(i)=a_noise(amp)
        enddo
        do i=1,N
          y(i)=y(i)+v*sin(angle(i))
          x(i)=x(i)+v*cos(angle(i))
          write(1,*) x(i),y(i),cos(angle(i)),sin(angle(i))
        enddo
        close(1)
      end

      subroutine initialize(N,L,PI,x,y,angle)
c       Initialize positions and orientations giving random values
        integer*4 N,L
        double precision PI
        double precision angle(N),x(N),y(N)
        double precision dran_u
        do i=1,N
          angle(i)=-PI+dran_u()*2.0d0*PI
          x(i)=dble(L)*dran_u()
          y(i)=dble(L)*dran_u()
        enddo
        return
      end

      double precision function a_noise(a)
c       Random number in [-a/2,a/2]
        double precision a
        double precision dran_u
        a_noise=-a/2.0d0+dran_u()*a
        return
      end
