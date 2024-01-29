!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!PARÁMETROS DE LA SIGMOIDE PARA EL VALLE (MESETA) DE INTENSIDAD!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ri_max = 5e-6; !en metros
ri_min = 5e-6; !en metros
z0i = 3.5e-3; !en metros
ki = 5e6; !en 1/metros
Li = ri_max-ri_min !PARÁMETRO AUXILIAR
ri(:) = ri_min + Li/(1+exp(ki*(this%spacez(:)-z0i))) !CALCULO DEL VECTOR SIGMOIDE

!PARÁMETROS DEL CANAL (100 micrometros de ancho del dominio y this%Nx = this%Ny = 100, con particiones dx = 1 micrometro).
Rc = 38e-6; !Radio del canal en metros
R0 = 80e-6; !Parámetro del canal en metros
Rv = 90e-6; !Parámetro del canal en metros
zshift = 0.75; !Posicion relativa del segundo pico de sobreionización (0 -> inicio del plasma, 0.5 -> mitad del canal, 1.0 -> final del canal)
cen_fac = 0.3; !Aumento de densidad electrónica sobre el valor sin sobreionizacion (cen_fac = 0.0) del segundo pico de sobreionización
sig_z = 0.25e-3; !Desviación estandar en dirección longitudinal (propagación) de la sobreionización del segundo pico (en metros)
sig_r = 15e-6; !Desviación estandar en dirección transversal (radial) de la sobreionización del segundo pico (en metros)
rlion_cen(:) = ri(:); !ALMACENA LA SIGMOIDE EN EL VECTOR RLION_CEN (antes constante)
z0_fac = 0.15; !Aumento de densidad electronica sobre el valor sin sobreionizacion (z0_fac = 0.0) al inicio de la columna de plasma
sig_z0 = 2e-3; !Desviación estandar en dirección longitudinal (propagación) de la sobreionización del primer pico (en metros)
sig_r0 = 15e-6; !Desviación estandar en dirección transversal (radial) de la sobreionización del primer pico (en metros)
rlion_z0 = 40.e-6; !(en metros)
          
z0_flag = 0.0;
cen_flag = 0.0;
          
if(z0_fac .gt. 0.0)then
  z0_flag = 1.0;
endif

if(cen_fac .gt. 0.0)then
  cen_flag = 1.0;
endif

do k = 1,this%Nx !Número de celdas en x
  do l = 1,this%Ny !Número de celdas en y
!INICIALIZACIÓN DEL CANAL DE PLASMA
#if 0
  !FORMA LONGITUDINAL DEL CANAL 
  if(this%pflag)then
    sigmaPz = this%plength/(2*(2*log(2.0))**(1/(2*this%nhgauss)))
    !N2 PERFIL DE DENSIDAD DE NEUTROS
    this%NZL3D(k,l,:) = exp(-0.5*((this%spacez(:)-this%pcenter)/sigmaPz)**(2*this%nhgauss))* &
      exp(-0.5*((k-this%Nx/2.0)/(this%Nx/5))**2-0.5*((l-this%Ny/2.0)/(this%Nx/5))**2)
  else
    !N2 PERFIL DE DENSIDAD DE NEUTROS
    this%NZL3D(k,l,:) = exp(-0.5*((k-this%Nx/2.0)/(this%Nx/2))**2-0.5* &
      ((l-this%Ny/2.0)/(this%Nx/2))**2)
  endif
    this%neprof3D(k,l,:) = this%NZL3D(k,l,:)
#endif           

!CANAL DE PLASMA INICIAL (asumo 100 micrometros de ancho del dominio).
#if 1
  !PREFORMA DEL CANAL DE PLASMA 
  radius = sqrt((this%spacex(k)-0.5*this%Lx)**2+(this%spacey(l)-0.5*this%Ly)**2)
  if(radius <= Rc)then
    this%NZL3D(k,l,:) = 1 + (radius/R0)**2
  else
    this%NZL3D(k,l,:) = (1+(Rc/R0)**2)*((Rv-radius)/(Rv-Rc))
  endif
      
  !DENSIDAD DE ELECTRONES 
  this%neprof3D(k,l,:) = this%NZL3D(k,l,:)*(1.0 + cen_fac*exp(-0.5*((this%spacez(:)-zshift*this%Lz)/sig_z)**2)* &
  exp(-0.5*((radius/sig_r)**2)) + z0_fac*exp(-0.5*(this%spacez(:)/sig_z0)**2)*exp(-0.5*(radius/sig_r0)**2))     
         
  !INTRODUCCIÓN DE LA SIGMOIDE DEL VALLE DE INTENSIDAD
  expion_cen(:) = exp(-0.5*(max(radius,rlion_cen(:))/sig_r)**2)
  expion_cen(:) = expion_cen(:)/exp(-0.5*(rlion_cen(:)/sig_r)**2)
  expion_z0 = exp(-0.5*(max(radius,rlion_z0)/sig_r0)**2)
  expion_z0 = expion_z0/exp(-0.5*(rlion_z0/sig_r0)**2)

  !DENSIDAD DE IONES DE Kr8+
  this%NZL3D(k,l,:) = this%NZL3D(k,l,:)*(1.0 - cen_flag*exp(-0.5*((this%spacez(:)-zshift*this%Lz)/sig_z)**2)* &
  exp(-0.5*(radius/sig_r)**2))*expion_cen*(1.0 - z0_flag*exp(-0.5*(this%spacez(:)/sig_z0)**2)*expion_z0) 
#endif
