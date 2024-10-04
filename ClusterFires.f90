! f2py3  -m ClusterFires -c ClusterFires.f90
!
!
subroutine cluster(x,y,z, s, frp, ipixel, icluster, idup, dupdist, N) 
    !
    ! Aggregates granule pixels into clusters
    
    implicit none

    real (kind=4), intent(in)    :: x(0:N-1), y(0:N-1), z(0:N-1), s(0:N-1), frp(0:N-1) 
                          ! Cartesian coordinates and sizes of pixels, pixelfrp
    integer (kind=8), intent(in) :: ipixel(0:N-1)! Pixel ID to mark duplicates
    integer (kind=4), intent(inout) :: icluster(0:N-1), idup(0:N-1) ! Cluster index, index of dup covered by the pixel
    real (kind=4), intent(inout) :: dupdist(0:N-1) !! Distance from duplicate
    integer (kind=4), intent(in) :: N
    !!Could not pass number of clusters..

    integer :: i,j, iold, inew, nclusters
    real (kind=4) :: dist2, s2
    REAL(kind=4), PARAMETER  :: F_NAN = TRANSFER(2143289344,1.0)


    icluster(0:N-1) = -1
    idup(0:N-1) = -1
    dupdist(0:N-1) = F_NAN
    nclusters=0

    do i = 0,N-1
      if (icluster(i) < 0) then
                icluster(i) = nclusters
                nclusters = nclusters + 1
      endif
      do j = i+1,N-1
          dist2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 
          s2 = (0.5*(s(i)+s(j)))**2
          if (dist2 < 2*s2 ) then ! adjacent?
             if (dist2 < 0.3*s2 ) then ! ## same fire, mark smaller FRP as duplicate
                    if (frp(i) > frp(j)) then
                      idup(j) = int(ipixel(i),kind=4)
                      dupdist(j) = sqrt(dist2)
                    else
                      idup(i) = int(ipixel(j),kind=4)
                      dupdist(i) = sqrt(dist2)
                    endif
             endif
             continue !! Fixme Disable clustering
             if (icluster(j) < 0) then !: ## Not yet assigned, just merge
               icluster(j) = icluster(i)
             else  !!! Reassign previous cluster
               iold = icluster(j)
               inew = icluster(i)
               where (icluster == iold) icluster = inew
             endif
           endif
      enddo
    enddo

end subroutine cluster


subroutine clustersOfFire(x,y,z,s, icluster, N, Nfires) 
    !
    ! Finds clusters that belong to the same fire
    ! HUOM! Sorted Z assumed, 
    !
    implicit none

    real (kind=4), intent(in)    :: x(0:N-1), y(0:N-1), z(0:N-1), s(0:N-1) ! Cartesian coordinates and size
    integer (kind=4), intent(inout) :: icluster(0:N-1) ! Cluster index
    integer (kind=4), intent(inout)  :: Nfires(1) !!!Number of aggregates
    integer (kind=4), intent(in) :: N

    integer :: i,j, iold, inew, nclusters
    real (kind=4) :: dist2, s2


    icluster(0:N-1) = -1
    nclusters=0

    do i = 0,N-1
      if (icluster(i) < 0) then
                icluster(i) = nclusters
                nclusters = nclusters + 1
      endif
      do j = i+1,N-1
          if (abs(z(i)-z(j)) > 30) exit !! 100 km should be far enough
          dist2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 
          s2 = (0.5*(s(i)+s(j)) )**2  
          if (dist2 < s2 ) then ! Same fire
             if (icluster(j) < 0) then !: ## Not yet assigned, just merge
               icluster(j) = icluster(i)
             else  !!! Reassign previous cluster
               iold = icluster(j)
               inew = icluster(i)
               where (icluster == iold) icluster = inew
             endif
          endif
      enddo
    enddo
    Nfires = nclusters

end subroutine clustersOfFire


 !*******************************************************

subroutine cluster2fire(tstamps, frps, hourvar, frpmean, N, asSILAM)
    implicit none
    !
    ! Fit dailymean FRP for observed values, summing corresponding FRP if timestamps coincide
    !  (likely from different clusters of the same satellite)
    ! 
    !  Algorithm approximately one from older SILAM fire source
    integer (kind=4), intent(in) :: tstamps(0:N-1)! Seconds since local daystart, SORTED!!!
    real (kind=4), intent(in)    :: frps(0:N-1), hourvar(0:23) ! Corresponding bserved FRP, prescribed variation
    real (kind=4), intent(inout)  :: frpmean(2) !! Effective FRP assuming landuse/falt diurnal var
    integer (kind=4), intent(in) :: N
    logical, intent (in) :: asSILAM

   
    real, parameter :: sigma_FRP_abs = 1. !MW ! 10 MW is a Blue-sky esimates by MAS (silam_v5_2, 18.11.2012) 
    real, parameter :: sigma_FRP_rel = 0.2 ! plus 20%
    real, parameter :: maxobservedfrac = 2. !Limit maximum FRP to 2xmax_observed 
    character(len = *), parameter :: sub_name = 'cluster2fire'


    integer :: i, j
    real :: sxx, sxy, sxxFlat, sxyFlat, FRP, rpcsig2

    sxx = 0.
    sxy = 0.
    sxxFlat = 0.
    sxyFlat = 0.
    FRP = 0.
    do i = 0, N-1
        FRP = FRP + frps(i)
        if (i < N-1 ) then
          if (tstamps(i) == tstamps(i+1)) then !! Another part of same fire,
                              !assume no same-timestamp detections from different satellites
            cycle
          endif
        endif
        j = tstamps(i) / 3600 !!Integer!!!
        rpcsig2 = 1./(sigma_FRP_abs*sigma_FRP_abs + FRP*FRP*sigma_FRP_rel*sigma_FRP_rel) ! 1/sigma2

        if (asSILAM) then
          ! As in SILAM
           sxx = sxx +  rpcsig2
           sxy = FRP*rpcsig2 /  hourvar(j)
         else
          !!!! To minimize  RMSE
           sxx = sxx + hourvar(j)*hourvar(j)*rpcsig2
           sxy = sxy + hourvar(j)*FRP*rpcsig2
         endif


        sxxFlat = sxxFlat + 1.
        sxyFlat = sxyFlat + FRP
        FRP = 0.
    enddo

    frpmean(1) = min(2*maxval(frps), sxy/sxx) !! Can't be twice the observed max
    frpmean(2) = sxyFlat/sxxFlat  !! Just mean obs FRP

end subroutine cluster2fire

 !*******************************************************

subroutine cluster2fireSILAM(tstamps, frps, hourvar, frpmean, N)
    implicit none
    !
    ! Fit dailymean FRP for observed values, summing corresponding FRP if timestamps coincide
    !  (likely from different clusters of the same satellite)
    ! 
    !  Algorithm approximately one from older SILAM fire source
    integer (kind=4), intent(in) :: tstamps(0:N-1)! Seconds since local daystart, SORTED!!!
    real (kind=4), intent(in)    :: frps(0:N-1), hourvar(0:23) ! Corresponding bserved FRP, prescribed variation
    real (kind=4), intent(inout)  :: frpmean(2) !! Effective FRP assuming landuse/falt diurnal var
    integer (kind=4), intent(in) :: N

    call  cluster2fire(tstamps, frps, hourvar, frpmean, N, .TRUE.)

end subroutine cluster2fireSILAM
 !*******************************************************

subroutine cluster2fireRMSE(tstamps, frps, hourvar, frpmean, N)
    implicit none
    !
    ! Fit dailymean FRP for observed values, summing corresponding FRP if timestamps coincide
    !  (likely from different clusters of the same satellite)
    ! 
    !  Algorithm approximately one from older SILAM fire source
    integer (kind=4), intent(in) :: tstamps(0:N-1)! Seconds since local daystart, SORTED!!!
    real (kind=4), intent(in)    :: frps(0:N-1), hourvar(0:23) ! Corresponding bserved FRP, prescribed variation
    real (kind=4), intent(inout)  :: frpmean(2) !! Effective FRP assuming landuse/falt diurnal var
    integer (kind=4), intent(in) :: N

    call  cluster2fire(tstamps, frps, hourvar, frpmean, N, .FALSE.)

end subroutine cluster2fireRMSE



