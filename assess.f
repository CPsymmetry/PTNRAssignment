	program neutron_populations
	implicit none
	character*30 files(3), uts
	character*30, allocatable :: isos(:)
	integer isonum, grmax, tm, iter
	integer i, j, k
	integer, allocatable :: intr(:)
	real tc, dndt, h
	real properties(3)
	real, allocatable :: gr(:,:,:), conc(:,:,:), npop(:,:)
	real, allocatable :: dcdt(:), beta(:), p(:), t(:)
c!	files: 1:proprties file, 2:input file, 3:output file
c!	uts: units for reactivity change. Nile and Dollar supported
c!	isos: isotope labels
c!	isonum: # of isotopes
c!	grmax: max # of groups for a single isotope
c!	tm: max time
c!	i,j,k: dummy variables for for loops
c!	intr: # of groups for each isotope
c!	tc: time constant - properties(1)
c!	dndt: differential of neutron population over time
c!	h: step size or tm/iter
c!	properties: contains all properties for input:
c!	1:time constant 2:initial reactivity 3:reactivitiy change
c!	gr: array of delayed neutron group constans for various fissile 
c!	    isotopes ordered by gr(isotope, group, decay constant (1) or 
c!	    Fractional Yield (2))
c!	conc: array of concentrations at time t
c!	npop: neutron population at index*step
c!	dcdt: differential of concentrations over time for each group
c!	beta: Total Yield for each isotope	
c!	p: reactivity
c!	t: time elapsed
	
	files(1) = 'properties.txt'
	
c!	retrieve properties defined in properties file.
	call get_properties(files,isonum,grmax,tm,uts,properties,iter)
	h = real(tm)/real(iter) !!step size
c!	allocate memory for arrays now that isonum and grmax have been retrieved
	allocate(gr(isonum, grmax, 2))
	allocate(isos(isonum))
	allocate(conc(isonum,grmax,iter))
	allocate(npop(isonum,iter))
	allocate(beta(isonum))
	allocate(intr(isonum))
	allocate(p(isonum))
	allocate(dcdt(grmax))
	allocate(t(iter))
c!	Get data for gr, isos, intr, and beta.
	call get_data(files(2), isonum, grmax, gr, isos, intr, beta)
	
c!	MAIN PROGRAM
	
c!	initial reactivity depends on what units. Could add a table of prefixes for milli, centi, etc.
	if (uts.eq.'Dollar') then
		p = properties(2)*beta
	elseif (uts.eq.'Nile') then
		do i = 1, isonum
			p(i) = properties(2)
		enddo
	endif
c!	calculating initial concentrations and setting neutron population to 1
	tc = properties(1) !!time constant for shorthand

	do i = 1, isonum
		npop(i,1) = 1
		t(1) = h
		do j = 1, intr(i)
			conc(i,j,1) = gr(i,j,2)/(tc*gr(i,j,1))
		enddo
	enddo			
c!	Loop over all isotopes and calculate values of concentration and
c!	neutron population over range of iter
	do i = 1, isonum
		do j = 2, iter
			t(j)  = real(j)*h
			
			do k = 1, intr(i)
				dcdt(k) = gr(i,k,2)*(npop(i,j-1)/tc) -
     +			gr(i,k,1)*conc(i,k,j-1)
     				conc(i,k,j) = conc(i,k,j-1) + h*dcdt(k)
     				if (conc(i,k,j).lt.0) then
     					conc(i,k,j) = 0.0
     				endif
      		enddo
	
			dndt = (p(i)-beta(i))*(npop(i,j-1)/tc) + 
     +					sum(gr(i,:,1)*conc(i,:,j-1))
     
			npop(i,j) = npop(i,j-1) + h*dndt
			
			if (npop(i,j).lt.0) then
				npop(i,j) = 0
			endif
				
c!			Calculate dc/dt and dn/dt at previous time. Step at
c!			small intervals then use eulers approximation for next
c!			step. This gives concentrations and neutron populations 
c!			at each interval of time.
c!			Eulers Approximation (1st order)
c!			y(n) = y(n-1) + h*y'(n-1)
		enddo
	enddo
	
	call output(files(3),npop,conc,t,isos,isonum,grmax,iter)
	
	end
	
	subroutine get_properties(files,isonum,grmax,tm,uts,properties,
     +									    	iter)
	character*30 files(3), uts
	real properties(3)
	integer isonum, grmax, tm
c!	gets properties isonum, grmax which is then used to allocate memory to arrays
	open(11, file=files(1), status='old')
	read(11, '(12X,A)') files(2)			!!input file
	read(11, '(18X,A)') files(3)			!!output file
	read(11, '(14X,I1)') isonum			!!num of isotopes
	read(11, '(12X,I2)') grmax			!!max num of groups
	read(11, '(10X,I3)') tm				!!max time		(seconds)
	read(11, '(12X,I4)') iter			!!number of iterations
	read(11, '(15X,F6.3)') properties(1)	!!time constant 	(seconds)	
	read(11, '(20X,F6.3)') properties(2)	!!initial reactivity
	read(11, '(19X,F6.3)') properties(3)	!!reactivity change
	read(11, '(18X,A)') uts				!!reactivity units
	close(11)
	end
	
	subroutine get_data(input, isonum, grmax, gr, isos, intr, beta)
	integer i, j, jx, group, isonum, grmax, intr(isonum)
	character*30  input, isos(isonum), iso
	real gr(isonum,grmax,2), beta(isonum), decay, yield	
c!	Gets input data and chucks it all into gr and isos arrays.
	
	open(21, file=input, status='old')
	
	do i = 1, isonum
		read(21,'(A5,1X,I1)') iso, jx
		isos(i) = iso
c!		print *, iso
		intr(i) = jx
		do j = 1, jx
			read(21,'(I1,1X,F6.4,1X,F9.7)') group, decay, yield
			gr(i,group,1) = decay
			gr(i,group,2) = yield
c!			print '(I1,1x,F6.4,1x,F9.7)', group, decay, yield
		enddo
	enddo
	
	do i = 1, isonum
		beta(i) = sum(gr(1,:,2))
	enddo
	
	close(21)
	
	end
	
	subroutine output(outfile,npop,conc,t,isos,isonum,grmax,iter)
	character*10 s1, isr
	logical bool
	integer isonum, grmax, iter
	integer i,j,k,n
	character*30 outfile, ftype, isos(isonum)
	real npop(isonum,iter), conc(isonum,grmax,iter), t(iter)
c!	collect all output data together into a readable fortmat and then
c!	save to the output file.
	
	ftype = '.txt'
	n = 0
	bool = .False.
	write(isr,'(I1)') isonum
	
c!	If the file name already exists, create a new filename with (n) on the end
	inquire(file=trim(outfile)//ftype,exist=bool)
	do while (bool)
		n = n + 1
		write(s1, '(I1)') n
		s1 = '(' // trim(s1) // ')'
		inquire(file=trim(outfile)//trim(s1)//ftype,exist=bool)
	enddo
c!	Concatenate the filename ready to open the file
	if (n.eq.0) then
		outfile = trim(outfile)//ftype
	else
		outfile = trim(outfile)//trim(s1)//ftype
	endif
	
	open(31,file=outfile, status='new')
c	write names of isotopes out
	write(31, '((A5,14X)'//isr//'(3X,A5,14X))') 
     +'time',(Isos(i), i=1,isonum)
c	write all data into corresponding collumns
	do j = 1, iter
		write(31, '((F7.4,12X),'//isr//'(F10.5,12X))')
     +	t(j),	(npop(i,j), i=1,isonum)
	enddo
	close(31)
	end
