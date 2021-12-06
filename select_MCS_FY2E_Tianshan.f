      program select_MSC_FY2E
      parameter(mis=-999,inum=1,nrecl=1201)
! Manually selection area
      parameter(yr1=2013,mn1=6,dy1=1,hr1=0) !Start time
      parameter(nyr=1,nmn=1,ndy=1,nhr=24,nmt=1) !Time counts
      parameter(icth=250) !TBB threshold
!---
      character*4 chyr(nyr)
      character*2 chmn(nmn),chdy(ndy),chhr(nhr),chmt(nmt)
      real lat(nrecl),lon(nrecl)
      logical alive1
      character ch(nrecl)*1
      real TBB(nmt,nrecl,nrecl)
      real,allocatable::TBB_flag(:,:,:),lat_flag(:),lon_flag(:)
      character*100 fin,fout
      
      real ns_ind(nmt,nrecl,nrecl),nsc(nmt,nrecl,nrecl)
      integer ns_flag(nmt,nrecl,nrecl)

!!!!  Note: for ifort, use '-assume byterecl' to compile this code
cccccc  TIME
ccc year
      do i=1,nyr
         ii=i+yr1-1
         write(chyr(i),"(i4)") ii
      end do
ccc month
      do i=1,nmn
         ii=i+mn1-1
         if(ii.lt.10) then
             write(chmn(i),"(i1)") ii
             chmn(i)="0"//trim(chmn(i))
         else
             write(chmn(i),"(i2)") ii
         end if
      end do
ccc day
      do i=1,ndy
         ii=i+dy1-1
         if(ii.lt.10) then
             write(chdy(i),"(i1)") ii
             chdy(i)="0"//trim(chdy(i))
         else
             write(chdy(i),"(i2)") ii
         end if
      end do
ccc hour
      do i=1,nhr
         ii=i+hr1-1
         if(ii.lt.10) then
             write(chhr(i),"(i1)") ii
             chhr(i)="0"//trim(chhr(i))
         else
             write(chhr(i),"(i2)") ii
         end if
      end do
ccc minute
      do i=1,nmt
         ii=(i-1)*30
         if(ii.lt.10) then
             write(chmt(i),"(i1)") ii
             chmt(i)="0"//trim(chmt(i))
         else
             write(chmt(i),"(i2)") ii
         end if
      end do
      NT=nmt
!---

cccccc  AREA
      lon_west=45.0          !!! Original FY-2E Lon-Lat area
      lon_east=165.0         !!! Lon: 45E ~ 165E
      lat_south=-60.0        !!! Lat: 60N ~ -60S
      lat_north=60.0
      grid_d=0.1                  !!! resolution 0.1x0.1 degree
      do i=1,nrecl
         lat(i)=lat_north-(i-1)*grid_d
         lon(i)=lon_west+(i-1)*grid_d
      end do

      latS=40.0              !!! Tianshan Area
      latN=50.0              !!! Lon: 76E ~ 96E
      lonW=76.0              !!! Lat: 40N ~ 50N
      lonE=96.0
      
      iy1=(lat_north-latS)/grid_d+1
      iy2=(lat_north-latN)/grid_d+1
      ix1=(lonW-lon_west)/grid_d+1
      ix2=(lonE-lon_west)/grid_d+1
      nx=ix2-ix1+1
      ny=iy1-iy2+1

      ALLOCATE(TBB_flag(nmt,ny,nx))
      ALLOCATE(lat_flag(ny),lon_flag(nx))
!---

cccccc  READING FILES
      do iyr=1,nyr
        do imn=1,nmn
          do idy=1,ndy
            do ihr=1,nhr
              TBB(:,:,:)=mis
              TBB_flag(:,:,:)=mis
              do imt=1,nmt
                fin='./DATA/FY2E_TBB_IR1_OTG_'//
!                fin='/data2/SATE/FY2E/FY2E_TBB_IR1/'//
!     &          chyr(iyr)//chmn(imn)//'/FY2E_TBB_IR1_OTG_'//
     &          chyr(iyr)//chmn(imn)//chdy(idy)//'_'//
     &          chhr(ihr)//chmt(imt)//'.AWX'
                print*, fin
                inquire(file=fin,exist=alive1)
                if(alive1) then
                    open(21,file='scale_Ic.txt')
                    open(22,file='index_Ic.txt')

                    open(11,file=fin,form='unformatted',
     &              status='old',access='direct',recl=nrecl)
                    do i=1,nrecl
                      read(11,rec=i+2) ch
                      do j=1,nrecl
                        TBB(imt,i,j)=real(ichar(ch(j))+100)
                        if(TBB(imt,i,j).le.100) TBB(imt,i,j)=mis
                       end do
                     end do
                     close(11)
                     print*, "min TBB=",minval(TBB,MASK=(TBB.ne.mis))
                     do iy=iy1,iy2,-1
                        iyy=iy1-iy+1
                       do ix=ix1,ix2
                         ixx=ix-ix1+1
                         TBB_flag(imt,iyy,ixx)=TBB(imt,iy,ix)
                         lat_flag(iyy)=lat(iy)
                         lon_flag(ixx)=lon(ix)
                       end do
                     end do

                     do iy=1,ny
                       do ix=1,nx
                        ! convection occurs when TBB<250K
                         if(TBB_flag(imt,iy,ix).lt.icth .and.
     &                   TBB_flag(imt,iy,ix).ne.mis) then
                            write(21,100) lon_flag(ix),lat_flag(iy)
                            write(22,110) ix,iy
                         end if
                       end do
                     end do
                     close(21)
                     close(22)

                     if(minval(TBB_flag(imt,:,:),
     &               MASK=(TBB_flag(imt,:,:).ne.mis)).lt.icth) then
!---

cccccc  CLUSTER
!using dbscan to cluster convective grids (TBB<250K) into one or several MCS
                        call system ('R CMD BATCH cluster_MCS.R')
                        open(31,file='index_Ic.txt')
                        open(32,file='data_temp_Ic.txt') !cluster result

                        ns_ind(imt,:,:)=mis
                        sk=0
                        n=0
                        do is=1,nrecl*nrecl
                          read(31,*,end=98)indx,indy
                          read(32,*,end=99)no_ind
                          ns_ind(imt,indy,indx)=no_ind
                          sk=sk+1
                        end do
98     continue
99     continue
                        close(31)
                        close(32)
                        n=int(maxval(ns_ind(imt,:,:)))
                        print*, 'Number of MCS groups:',n
                        do ntt=1,n
                          mm=0
                          do ix=1,nx
                              do iy=1,ny
                                  if(ns_ind(imt,iy,ix).eq.ntt
     &                            .and.TBB_flag(imt,iy,ix).lt.icth .and.
     &                            TBB_flag(imt,iy,ix).ne.mis)
     &                            mm=mm+1
                              end do
                          end do

                          do ixx=1,nx
                            do iyy=1,ny
                              if(ns_ind(imt,iyy,ixx).eq.ntt
     &                        .and.TBB_flag(imt,iy,ix).lt.icth .and.
     &                        TBB_flag(imt,iy,ix).ne.mis)
     &                        nsc(imt,iyy,ixx)=mm

                              if(nsc(imt,iyy,ixx).gt.40) then
                                  ns_flag(imt,iyy,ixx)=
     &                            ns_ind(imt,iyy,ixx)
                              else
                                  ns_flag(imt,iyy,ixx)=0
                              end if
                            end do
                          end do

                        end do
                     end if

                else !fin not exist
                    TBB(imt,:,:)=mis
                    write(*,*) trim(fin)//" not exist!"
                end if
              end do    !end of minute loop

              if(any(ns_flag.gt.0)) then
                open(41,file='temp_info.txt')
                write(41,*)chyr(iyr),chmn(imn),chdy(idy),
     &          chhr(ihr),chmt(1)
                write(41,*)lon_flag(1),lon_flag(nx),
     &          lat_flag(1),lat_flag(ny)
                write(41,*)nx,ny,NT
                close(41)
                open(42,file='temp.txt')

                do imt=1,nmt
                  write(42,*)((ns_flag(imt,iy,ix),ix=1,nx),iy=1,ny)
                  print*,minval(nsc,MASK=(nsc.ne.mis)),
     &            maxval(nsc,MASK=(nsc.ne.mis))
                  ns_flag(imt,:,:)=mis
                  nsc(imt,:,:)=mis
                end do
                close(42)
!---

cccccc  OUTPUT
!write into NetCDF data & draw the pic
                call system ('ncl write_MCS_to_ncdata.ncl')
!---
              end if

            end do    !end of hour loop
          end do    !end of day loop
        end do    !end of month loop
      end do    !end of year loop
      

cccccc  FINALIZE
      deallocate(TBB_flag,lat_flag,lon_flag)
      call system ('rm -f ./*txt')
!---
      
100   format(2f11.3)
110   format(2i7)

      
      end
