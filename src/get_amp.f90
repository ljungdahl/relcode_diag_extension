subroutine get_amp(pert_wave,size_rhs,Bbspline1_mat,Bbspline2_mat,dBbspline1_mat,dBbspline2_mat,pos_vec, &
pvs,size1,size2,omega,w_in,sw,eig_val_mat,l_vec,slv,amp,phaseF,phaseG,pcurF,pcurG,pcur,kappa_vec,skv,fsc_inv, &
z,energy,rpa_conv)
implicit none

integer,intent(in)::size_rhs,size1,size2,pvs,sw,slv,skv
integer,dimension(slv),intent(in)::l_vec
integer,dimension(skv),intent(in)::kappa_vec
integer,dimension(1,2,sw),intent(in)::w_in
complex*16,intent(in)::omega
complex*16,dimension(pvs),intent(in)::pos_vec
complex*16,dimension(pvs,size1+2),intent(in)::Bbspline1_mat,dBbspline1_mat
complex*16,dimension(pvs,size2+2),intent(in)::Bbspline2_mat,dBbspline2_mat
complex*16,dimension(size_rhs,3,sw),intent(in)::pert_wave
complex*16,dimension(size_rhs-size_rhs/2-1,skv),intent(in)::eig_val_mat
complex*16,dimension(5,3,sw),intent(out)::pcurF,pcurG,pcur
complex*16,dimension(sw),intent(out)::energy
real*8,dimension(5,3,sw),intent(out)::amp,phaseF,phaseG
real*8,intent(in)::fsc_inv,z
logical,dimension(sw,3),intent(in)::rpa_conv

integer,dimension(5)::pos_ind
complex*16::clmb1,clmb2,dclmb1,dclmb2,clmb_tot,sig,negi
complex*16,dimension(3)::rhoF,rhoG
complex*16,dimension(5,3)::tmp,dtmp
integer::k,n

!z = 1.d0
negi = (0.d0,1.d0)

amp 	= 0.d0
phaseF 	= 0.d0
phaseG 	= 0.d0
energy 	= (0.d0,0.d0)
pcur 	= (0.d0,0.d0)

k = 1
do while(real(pos_vec(k)) < 1.1d0/3.d0*real(pos_vec(size(pos_vec))))
	k = k+1
end do
pos_ind = (/k,k+75,k+150,k+225,k+300/)

!!!DEBUG!!!
!write(*,*)'Positions:',pos_vec(pos_ind)
!write(*,*)' '
!!!DEBUG END !!!

do k=1,sw
if (rpa_conv(k,3))then

	!!!! TEST OUTPUT !!!!
	!write(*,*)'l in get_amp ',l_vec(w_in(1,2,k))+1
	!write(*,*)'lr_charge in get_amp ',z
	!!!! END TEST OUTPUT !!!!

	energy(k) = omega+eig_val_mat(w_in(1,1,k),w_in(1,2,k))

	if (real(energy(k)) .ge. -15.d0)then
	do n=1,5

		rhoF = matmul(Bbspline1_mat(pos_ind(n),2:size1+1),pert_wave(1:size1,:,k))
		rhoG = matmul(Bbspline2_mat(pos_ind(n),2:size2+1),pert_wave(size1+1:size_rhs,:,k))

		if (l_vec(w_in(1,2,k)).ne.0)then
			call coul_gen_cc_rel(z,energy(k),real(pos_vec(pos_ind(n))),l_vec(w_in(1,2,k))-1,.false.,clmb1, &
					clmb2,dclmb1,dclmb2,sig)

			!!!! TESTING different clmb_tot setup !!
			clmb_tot = (clmb2+negi*clmb1)

			!clmb_tot = -negi*(clmb2+negi*clmb1)		!!! Temp taken out (might be more then temp

			!clmb_tot = (clmb1-negi*clmb2)
			!!!! TESTING
			!write(*,*)'l = ',l_vec(w_in(1,2,k))-1,' C_phase = ',sig
			!!!! TESTING DONE


			if (kappa_vec(w_in(1,2,k)) < 0)then
				amp(n,1,k) 		= (abs(rhoF(1))+abs(rhoG(1)))/abs(clmb_tot)
				phaseF(n,1,k) 	= atan2(aimag(rhoF(1)),real(rhoF(1)))-atan2(aimag(clmb_tot),real(clmb_tot))
				phaseG(n,1,k) 	= atan2(aimag(-negi*rhoG(1)),real(-negi*rhoG(1))) &
									-atan2(aimag(clmb_tot),real(clmb_tot))
			else
				amp(n,1:2,k) 	= (abs(rhoF(1:2))+abs(rhoG(1:2)))/abs(clmb_tot)
				phaseF(n,1:2,k) = atan2(aimag(rhoF(1:2)),real(rhoF(1:2))) &
									-atan2(aimag(clmb_tot),real(clmb_tot))
				phaseG(n,1:2,k) = atan2(aimag(-negi*rhoG(1:2)),real(-negi*rhoG(1:2))) &
									-atan2(aimag(clmb_tot),real(clmb_tot))
			end if
		end if

		! NOTE(anton): The appropriate Coulomb-function that is a part of the perturbed wave function is calculated here.
		! It is then used to get the non-Coulomb amplitude and the phase relative to the Coulomb phase, of the perturbed
		! wave function.
		call coul_gen_cc_rel(z,energy(k),real(pos_vec(pos_ind(n))),l_vec(w_in(1,2,k))+1,.false.,clmb1,clmb2, &
				dclmb1,dclmb2,sig)

		!!!! TESTING different clmb_tot setup !!
		clmb_tot = (clmb2+negi*clmb1)

		!clmb_tot = -negi*(clmb2+negi*clmb1)		!!! Temp taken out (might be more then temp


		!clmb_tot = (clmb2-negi*clmb1)
		!!!! TESTING
		!write(*,*)'l = ',l_vec(w_in(1,2,k))+1,' C_phase = ',sig
		!!!! TESTING DONE

		if (kappa_vec(w_in(1,2,k)) < 0)then
			amp(n,2:3,k) 	= (abs(rhoF(2:3))+abs(rhoG(2:3)))/abs(clmb_tot)
			phaseF(n,2:3,k) = atan2(aimag(rhoF(2:3)),real(rhoF(2:3)))-atan2(aimag(clmb_tot),real(clmb_tot))
			phaseG(n,2:3,k) = atan2(aimag(-negi*rhoG(2:3)),real(-negi*rhoG(2:3))) &
								-atan2(aimag(clmb_tot),real(clmb_tot))

		else
			amp(n,3,k) 		= (abs(rhoF(3))+abs(rhoG(3)))/abs(clmb_tot)
			phaseF(n,3,k) 	= atan2(aimag(rhoF(3)),real(rhoF(3)))-atan2(aimag(clmb_tot),real(clmb_tot))
			phaseG(n,3,k) 	= atan2(aimag(-negi*rhoG(3)),real(-negi*rhoG(3))) &
								-atan2(aimag(clmb_tot),real(clmb_tot))
		end if

		pcur(n,:,k)	= 2*negi*(conjg(rhoF)*rhoG-conjg(rhoG)*rhoF)
!!!!!!!DEBUG !!!!!!!!!!!!!
!		write(*,*)'k = ',k,' n = ',n
!		write(*,*)'pcur: ',pcur(n,:,k)
!		write(*,*)' '
!!!!!!!DEBUG END !!!!!!!!!!!!!

	end do
	end if

	tmp 			= matmul(Bbspline1_mat(pos_ind,2:size1+1),pert_wave(1:size1,:,k))
	dtmp 			= matmul(dBbspline1_mat(pos_ind,2:size1+1),pert_wave(1:size1,:,k))
	pcurF(:,:,k) 	= negi*(conjg(tmp)*dtmp-tmp*conjg(dtmp))
	tmp 			= matmul(Bbspline2_mat(pos_ind,2:size2+1),pert_wave(size1+1:size_rhs,:,k))
	dtmp 			= matmul(dBbspline2_mat(pos_ind,2:size2+1),pert_wave(size1+1:size_rhs,:,k))
	pcurG(:,:,k) 	= negi*(conjg(tmp)*dtmp-tmp*conjg(dtmp))

end if
end do

!!!!!! DEBUG !!!!!!!
!write(*,*)pcurF(:,:,1)
!write(*,*)pcur(:,:,1)
!!!!!! DEBUG END!!!!!!!

!write(*,*)'Position on grid for amplitude and phase of rho:'
!write(*,*)real(pos_vec(pos_ind))
!write(*,*)' '

!!!! TEST TEMP
!write(*,*)'Energy out get_amp ',energy
!!!! END TEST


end subroutine
