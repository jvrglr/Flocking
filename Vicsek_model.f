      program Vicksek_model
c        atan2 documentation--> https://gcc.gnu.org/onlinedocs/gcc-4.2.3/gfortran/ATAN2.html
c        Quick presentation of Link Cell Method -->https://en.wikipedia.org/wiki/Cell_lists
c        REFERENCES:
C       ---------------------------------------------------------------------------
C        -->Vicsek, T., Czirók, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995).
C        Novel type of phase transition in a system of self-driven particles.
C        Physical review letters, 75(6), 1226.
C       ---------------------------------------------------------------------------
        implicit none
        integer*4 N,L,t,ncel,total_cell,realiz,steps
        double precision PI,mean,s_m,c_m,r,v,d,rbird,h
        !d=length of cell for Link cell method.There will be (L/a)**2 cells
        !r=interaction length
        parameter (N=1000,L=20,realiz=500
     c   ,PI=4.0d0*datan(1.0d0),r=1.0d0,v=1.0d-1,d=r) !Number of birds
        parameter (ncel=int(L/d),total_cell=ncel**2+4*(ncel+1))
        double precision angle(2*N),new_angle(N) !orientation of each bird, angle in [-PI,PI]
        double precision x(2*N),y(2*N) !x,y position of each bird
        double precision a_noise,amp,amp_min,amp_max !angular noise and amplitude of noise
        double precision dran_u,dist,v_mean
        integer*4 cell(total_cell,N)
        integer*4 info_c(total_cell)
        integer*4 box(N)
        integer*4 neigh(ncel**2,8)
        integer*4 i,j,k,aux,bird,count
        double precision mag,c,m2,new,old,error,dev,tau
        real start,finish !To compute CPU time
        character(len=20) str

c       JUST RUN ONCE------------------------------------
        call dran_ini(1994) !Seed for random number generator
        call initialize(N,L,PI,x,y,angle) ! Initial condition
        call create_neigh(ncel,neigh) ! Neighbors of each cell
        call cpu_time(start)
c       -------------------------------------------
        !Amp is not a parameter since it's value is changed in the program
        !when measuring phase diagram
        amp=0.4 !Noise amplitude

        steps=40 !number of points
        amp_max=5.0
        amp_min=0.0
        h=(amp_max-amp_min)/steps
          t=0
          do t=1,300
            open(unit=1, file='picture.'//trim(str(t))//'.dat',
     c            status="unknown")
           do i=1,N
             write(1,*) x(i),y(i),2*v*cos(angle(i)),
     c    2*v*sin(angle(i))
           enddo
            !In which cell is each bird?+periodic boundary conditions
            call update_cells(
     c        N,L,d,ncel,box,x,y,info_c,total_cell,cell)

            !Create/update ghost cells
            call update_ghost(
     c        x,y,info_c,ncel,N,total_cell,L,cell,angle)

            !Update aligment
            call interaction(
     c  N,total_cell,ncel,angle,new_angle,neigh,info_c,x,y,box,r,cell)

              !UPDATE POSITIONS
              do i=1,N
                angle(i)=new_angle(i)+a_noise(amp)
                y(i)=y(i)+v*sin(angle(i))
                x(i)=x(i)+v*cos(angle(i))
              enddo
            close(1)
          enddo

        call cpu_time(finish)
        write(*,*) 'Time in seconds',finish-start


      end

      subroutine interaction(
     c  N,total_cell,ncel,angle,new_angle,neigh,info_c,x,y,box,r,cell)
       integer*4 N,L,t,ncel,total_cell
       double precision PI,mean,s_m,c_m,r,v,d,rbird
       double precision angle(2*N),new_angle(N) !orientation of each bird, angle in [-PI,PI]
       double precision x(2*N),y(2*N) !x,y position of each bird
       double precision a_noise,amp !angular noise and amplitude of noise
       double precision dran_u,dist,v_mean
       integer*4 cell(total_cell,N)
       integer*4 info_c(total_cell)
       integer*4 box(N)
       integer*4 neigh(ncel**2,8)
       integer*4 i,j,k,aux,bird,count,realiz

       do i=1,N
         count=1 !numbers of birds within radious r
         s_m=sin(angle(i))
         c_m=cos(angle(i))
           do j=1,8
             aux=neigh(box(i),j) !label of adjacen cell
             do k=1,info_c(aux) !pass through all birds in adjacent cell
               bird=cell(aux,k)
               rbird=dist(x(i),x(bird),y(i),y(bird))
               if (rbird.le.r) then !if birds are close enough they interact
                 count=count+1
                 s_m=s_m+sin(angle(bird))
                 c_m=c_m+cos(angle(bird))
               endif
             enddo
           enddo
           s_m=s_m/dble(count)
           c_m=c_m/dble(count)
           new_angle(i)=atan2(s_m,c_m)
         enddo

          return
        end

      character(len=20) function str(k)
      !   "Convert an integer to string."
          integer, intent(in) :: k
          write (str, *) k !write to a string
          str = adjustl(str)
      end function str

      subroutine update_cells(
     c      N,L,d,ncel,box,x,y,info_c,total_cell,cell)
c       Assign cell to each bird
c       Implement Periodic boundary conditions
        integer*4 N,L,total_cell,ncel
        double precision d
        double precision x(2*N),y(2*N)
        integer*4 box(N),info_c(total_cell),cell(total_cell,N)
        integer*4 label,xcell,ycell,i
        info_c=0

c       Boundary conditions----------------
        do i=1,N
          if (x(i).ge.L) then
            x(i)=x(i)-L
          else if (x(i).lt.0) then
            x(i)=x(i)+L
          endif
          if (y(i).lt.0) then
            y(i)=y(i)+L
          else if (y(i).ge.L) then
            y(i)=y(i)-L
          endif
c         -------------------------------------
          !compute label of cell in which is bird x(i),y(i)
          xcell=int(x(i)/d)
          ycell=int(y(i)/d)
          label=ycell*ncel+xcell+1 !add +1 to have label in [1,ncel^2]
          box(i)=label
          info_c(label)=info_c(label)+1
          cell(label,info_c(label))=i

        enddo

        return
      end

      subroutine initialize(N,L,PI,x,y,angle)
c       Initialize positions and orientations randomly
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

      double precision function v_mean(N,angle)
        integer*4 N,i
        double precision angle(2*N)
        double precision sumx,sumy

        sumx=0.0d0
        sumy=0.0d0
        do i=1,N
          sumx=sumx+cos(angle(i))
          sumy=sumy+sin(angle(i))
        enddo
        v_mean=dsqrt(sumx**2.0d0+sumy**2.0d0)/dble(N)
        return
      end

      double precision function dist(x1,x2,y1,y2)
        !disance between points (x1,x2) and (y1,y2)
        !dist(x,0,y,0)=module(x,y)
        double precision x1,x2,y1,y2
        dist=dsqrt((x1-x2)**2.0d0+(y1-y2)**2.0d0)
        return
      end

      subroutine create_neigh(ncel,neigh)
c       >Store label of adjacent cells for each cell
c       >Periodic boundary conditions have been taken into account
c       >It could be computionally expensive, but is only
!       runned once
        integer*4 ncel,i,x,y
        integer*4 neigh(ncel**2,8)
        do i=1,ncel**2
          x=mod((i-1),ncel)
          y=int((i-1)/ncel)
          neigh(i,1)=x+1+y*ncel+1 !right neighbor (NN)
          neigh(i,2)=x-1+y*ncel+1 !left neighbor (NN)
          neigh(i,3)=x+1+(y+1)*ncel+1 !right+up neigbor (diagonal)
          neigh(i,4)=x-1+(y-1)*ncel+1 !left+down neigh (diagonal)
          neigh(i,5)=x+1+(y-1)*ncel+1 !right+down neigh (diagonal)
          neigh(i,6)=x-1+(y+1)*ncel+1 !left +up neigh (diagonal)
          neigh(i,7)=x+(y-1)*ncel+1 ! down neigh (NN)
          neigh(i,8)=x+(y+1)*ncel+1 ! up neigh (NN)


          !Boundary conditions: Labeling gohst cells on arista
          if (y.eq.(ncel-1)) then
            neigh(i,8)=ncel**2+x+1
            neigh(i,3)=ncel**2+x+2
            neigh(i,6)=ncel**2+x
          endif
          if (y.eq.0) then
            neigh(i,7)=ncel**2+ncel+x+1
            neigh(i,5)=ncel**2+ncel+x+2
            neigh(i,4)=ncel**2+ncel+x
          endif
          if (x.eq.(ncel-1)) then
            neigh(i,1)=ncel**2+2*ncel+y+1
            neigh(i,3)=ncel**2+2*ncel+y+2
            neigh(i,5)=ncel**2+2*ncel+y
          endif
          if (x.eq.0) then
            neigh(i,2)=ncel**2+3*ncel+y+1
            neigh(i,6)=ncel**2+3*ncel+y+2
            neigh(i,4)=ncel**2+3*ncel+y
          endif
        enddo

        !Boundary conditions: Labeling gohst cells on vertex
        neigh(ncel**2,3)=ncel**2+4*ncel+1
        neigh(ncel,5)=ncel**2+4*ncel+2
        neigh(1,4)=ncel**2+4*ncel+3
        neigh(ncel*(ncel-1)+1,6)=ncel**2+4*ncel+4
        return
      end

      subroutine update_ghost(
     c  x,y,info_c,ncel,N,total_cell,L,cell,angle)
c      >Create ghost cells as copy of cells at the boundaries
c      >Since we have to work with the labels of the cells, it
c      can be messy code. However it isreally fast.
       integer*4 N,ncel,total_cell,L,box,up,down,bird
       integer*4 info_c(total_cell),cell(total_cell,N)
       double precision x(2*N),y(2*N),angle(2*N)
       integer*4 i,aux
       box=ncel**2
       up=box!I create up to make code more readable
       down=ncel**2-ncel

       do i=1,ncel
          !Ghost cells above y=L
         aux=info_c(i)!copy first ncel cells, which corresponds to y=0
         info_c(up+i)=aux
         do j=1,aux
           bird=cell(i,j)
           cell(up+i,j)=bird+N
           x(bird+N)=x(bird)
           y(bird+N)=y(bird)+L
           angle(bird+N)=angle(bird)
         enddo
         !Ghost cells below y=0
         aux=info_c(down+i) !copying y=ncel-1 cells
         info_c(box+ncel+i)=aux
         do j=1, aux
           bird=cell(down+i,j)
           cell(box+ncel+i,j)=bird+N
           x(bird+N)=x(bird)
           y(bird+N)=y(bird)-L
           angle(bird+N)=angle(bird)
         enddo

         !Ghost cells at right of x=ncel-1
         aux=info_c((i-1)*ncel+1)!copying x=0 cells
         info_c(box+2*ncel+i)=aux
         do j=1, aux
           bird=cell((i-1)*ncel+1,j)
           cell(box+2*ncel+i,j)=bird+N
           x(bird+N)=x(bird)+L
           y(bird+N)=y(bird)
           angle(bird+N)=angle(bird)
         enddo

         !Ghost cells at left of x=0
         aux=info_c(i*ncel) !copy cell at x=L
         info_c(box+3*ncel+i)=aux
         do j=1, aux
           bird=cell(i*ncel,j)
           cell(box+3*ncel+i,j)=bird+N
           x(bird+N)=x(bird)-L
           y(bird+N)=y(bird)
           angle(bird+N)=angle(bird)
         enddo
       enddo

       !Ghost cell in x=ncel y=ncel
       aux=info_c(1) !Copy cell x=0 y=0
       info_c(ncel**2+4*ncel+1)=aux
       do j=1, aux
         bird=cell(1,j)
         cell(ncel**2+4*ncel+1,j)=bird+N
         x(bird+N)=x(bird)+L
         y(bird+N)=y(bird)+L
         angle(bird+N)=angle(bird)
       enddo

       !Ghost cell in x=ncell, y=-1
       aux=info_c(box-ncel+1) !copy cell x=0,y=ncell-1
       info_c(box+4*ncel+2)=aux
       do j=1, aux
         bird=cell(box-ncel+1,j)
         cell(box+4*ncel+2,j)=bird+N
         x(bird+N)=x(bird)+L
         y(bird+N)=y(bird)-L
         angle(bird+N)=angle(bird)
       enddo

       !Ghost cell in x=-1,y=-1
       aux=info_c(ncel**2) !copy cell in x=Y=ncel-1
       info_c(ncel**2+4*ncel+3)=aux
       do j=1, aux
         bird=cell(ncel**2,j)
         cell(box+4*ncel+3,j)=bird+N
         x(bird+N)=x(bird)-L
         y(bird+N)=y(bird)-L
         angle(bird+N)=angle(bird)
       enddo

       !Ghost cell in y=ncel x=-1
       aux=info_c(ncel) !copy cell in y=0,x=ncell-1
       info_c(total_cell)=aux
       do j=1, aux
         bird=cell(ncel,j)
         cell(total_cell,j)=bird+N
         x(bird+N)=x(bird)-L
         y(bird+N)=y(bird)+L
         angle(bird+N)=angle(bird)
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
