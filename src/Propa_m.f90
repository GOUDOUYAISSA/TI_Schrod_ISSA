module Propa_m
  USE NumParameters_m
  implicit none


contains
  SUBROUTINE propagation(psif,psi0,H,t0,tf,delta_t)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psif
  TYPE (psi_t),  intent(in)    :: psi0
  TYPE(Op_t),    intent(in)    :: H

  real(kind=Rk), intent(in)    :: t0,tf,delta_t

  ! variables locales
  real(kind=Rk) :: t,t_deltat
  integer       :: i,nt
  TYPE (psi_t)  :: psi,psi_dt

  write(out_unitp,*) 'BEGINNIG propagation', t0,tf,delta_t

  nt = int((tf-t0)/delta_t)

  CALL init_psi(psi,psi0%Basis,cplx=.TRUE.) ! to be changed
  CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.) ! to be changed

  psi%CVec(:) = psi0%CVec

  DO i=0,nt-1

    t = i*delta_t
    t_deltat = t + delta_t

    write(out_unitp,*) 'march taylor',i,t,t_deltat
    CALL march_taylor(psi,psi_dt,H,t,t_deltat)

    psi%CVec(:) = psi_dt%CVec

  END DO

  psif%CVec(:) = psi%CVec

  CALL dealloc_psi(psi)
  CALL dealloc_psi(psi_dt)


  write(out_unitp,*) 'END propagation'

  END SUBROUTINE propagation

  SUBROUTINE march_taylor(psi,psi_dt,H,t,t_deltat)
  USE op_m
  USE psi_m
 
  TYPE (psi_t),  intent(inout) :: psi_dt
  TYPE (psi_t),  intent(inout)    :: psi
  TYPE (psi_t)                 :: Hpsi
  TYPE(Op_t),    intent(in)    :: H

  real(kind=Rk), intent(in)    :: t,t_deltat
  real(kind=Rk)                :: delta_t
  ! variables locales
  real(kind=Rk)                :: Rkk,Norm 
  integer       :: kk
     

  write(out_unitp,*) 'BEGINNIG march_taylor',t,t_deltat
  write(out_unitp,*) 'psi',psi%CVec

  psi_dt%CVec    = psi%CVec ; Rkk =1
  !======================================ordre 1==========================
  !Rj = 1
    !CALL calc_OpPsi(H,psi,psi)
    !Rj = Rj*delta_t
    !psi_dt%CVec(:) = psi_dt%CVec(:) + -EYE*delta_t*psi%CVec(:)
   !======================================================================
   
                              !ordre 2 et plus
        Do kk = 1,400,1                      
            Rkk = Rkk*(delta_t/kk)
               CALL calc_OpPsi(H,psi,Hpsi)
                     Hpsi%CVec(:)   = -EYE*Hpsi%CVec(:)
                     psi_dt%CVec(:) = psi_dt%CVec(:) + Rkk*Hpsi%CVec(:) 
                        psi%CVec(:) = Hpsi%CVec(:)
                            Norm    = Rkk*sqrt(abs(dot_product(psi%CVec,psi%CVec)))
                            
                                  IF (Norm .lt.10**-10) then
                                       exit
                                          write(out_unitp,*) 'taylor condition is fulfild for kk=',kk
                                  Endif           
                          
        End do

  write(out_unitp,*) 'psi_dt',psi_dt%CVec
  write(out_unitp,*) 'END march_taylor'

  END SUBROUTINE march_taylor

end module Propa_m
